import re
import torch
from torch import nn
import numpy as np
import os
import subprocess
import torch
from transformers import XLNetTokenizer, XLNetLMHeadModel
from .protein_mpnn.protein_mpnn_utils import ProteinMPNN
import pickle as pkl

from pathlib import Path
from tr_rosetta_pytorch import trRosettaNetwork
from tr_rosetta_pytorch.utils import preprocess
from tr_rosetta_pytorch.cli import DEFAULT_MODEL_PATH

from .utils import AMINO_ACID_ORDER

class XLNetWrapper(nn.Module):
    def __init__(self, model_name='Rostlab/prot_xlnet', device=None):
        super().__init__()
        if device is not None:
            self.device = device
        elif torch.cuda.is_available():
            self.device = torch.device("cuda:0")
        else:
            self.device = torch.device('cpu')

        self.model = XLNetLMHeadModel.from_pretrained(model_name)
        self.model.to(self.device)
        self.model.eval()
        self.tokenizer = XLNetTokenizer.from_pretrained(model_name)

        xlnet_vocab_dict = self.tokenizer.get_vocab()
        xlnet_vocab_dict['▁X'] = xlnet_vocab_dict['X']
        self.canonical_idx_to_xlnet_idx = torch.tensor([xlnet_vocab_dict['▁' + aa] for aa in AMINO_ACID_ORDER])

    def forward(self, seq, decode_order, token_to_decode, struct=None, mask_type='bidirectional_autoregressive', temperature=1.0):
        """Accept an amino acid sequence, return class probabilities for the next token
        Args:
            seq (len N list of len L_seq str): a string representation of an amino 
                acid sequence with unknown residues indicated with a dash (-)
            decode_order (len L list): list of the order of indices to decode.
                This determines the values in the permutation mask. Each index 
                attends to all indices that occur previously in the decoding_order.
            token_to_decode (len N tensor): index in the range [0, L-1] indicating which 
                token to predict for each item in the batch.
        Returns:
            probs ((20) torch.Tensor): a vector of probabilities for the next
                token
        """        
        # Replace rare amino acids with "X"
        seq = [re.sub(r"[UZOB]", "X", s) for s in seq]
        # Huggingface XLNet expects a space-separated sequence
        seq = [" ".join(s) for s in seq]
        # Mask '-' tokens
        seq = [re.sub(r"-", "<mask>", s) for s in seq]
        input_ids = self.tokenizer(seq, return_tensors='pt', add_special_tokens=False)['input_ids']
        seq_len = input_ids.shape[-1]
        
        n_tokens = len(token_to_decode)
        
        # perm_mask should be the mask for query stream attention (don't allow items to see self), because the xlnet 
        # implementation subtracts an identity matrix from perm_mask to get the content stream attention, where items are
        # masked from seeing self. The target_mapping ensures that we use query stream attention when predicting
        # the target elements https://github.com/huggingface/transformers/blob/v4.23.1/src/transformers/models/xlnet/modeling_xlnet.py#L1172
        # In query stream attention you should always be masked from yourself because tokens are always masked from themselves during training
        perm_mask = torch.ones((n_tokens, seq_len, seq_len), dtype=torch.float) # perm_mask[0, j, k] = 1 means that the jth token cannot see the kth token
        if mask_type == 'unidirectional_autoregressive': # Allow each token to see tokens preceding it in decode order
            for i, tok in enumerate(token_to_decode):
                token_to_decode_idx = decode_order.index(tok)
                # +1 because we want to include the tokens that the token_to_decode can see
                for j, idx in enumerate(decode_order[:token_to_decode_idx+1]):
                    perm_mask[i, idx, decode_order[:j]] = 0.0
            # decode_order: [2, 1, 0]
            # i = 0 -> decode pos 2
            #   j = 0, idx = 2
            #           [1, 1, 1]
            #           [1, 1, 1]
            #           [1, 1, 1]
            # i = 1 -> decode pos 1
            #   j = 0, idx = 2
            #   j = 1, idx = 1
            #       perm_mask[1, 1, [2]] = 0.
            #           [1, 1, 1]
            #           [1, 1, 0]
            #           [1, 1, 1]
            # i = 2 -> decode pos 0
            #   j = 0, idx = 2
            #   j = 1, idx = 1
            #       perm_mask[2, 1, [2]] = 0.
            #           [1, 1, 1]
            #           [1, 1, 0]
            #           [1, 1, 1]
            #   j = 2, idx = 0
            #       perm_mask[2, 0, [2, 1]] = 0.
            #           [1, 0, 0]
            #           [1, 1, 0]
            #           [1, 1, 1]
        elif mask_type == 'bidirectional_autoregressive':
            for i, tok in enumerate(token_to_decode):
                token_to_decode_idx = decode_order.index(tok)
                # Iterate over all tokens up to and including tok
                for prev_tok_1 in decode_order[:token_to_decode_idx + 1]:
                    # Iterate over all tokens preceding tok
                    for prev_tok_2 in decode_order[:token_to_decode_idx]:
                        if prev_tok_1 == prev_tok_2:
                            # Never let a token see itself in query-stream attention
                            continue
                        # Allow all tokens up to and including tok to see tokens preceding tok
                        perm_mask[i, prev_tok_1, prev_tok_2] = 0.0

            #       [1, 1, 1, 1]
            #       [1, 1, 1, 1]
            #       [1, 1, 1, 1]
            #       [1, 1, 1, 1]
            #
            #       [1, 1, 1, 1]
            #       [1, 1, 1, 1]
            #       [1, 1, 1, 0]
            #       [1, 1, 1, 1]
            #
            #       [1, 1, 1, 1]
            #       [1, 1, 0, 0]
            #       [1, 1, 1, 0]
            #       [1, 1, 0, 1]
            #
            #       [1, 0, 0, 0]
            #       [1, 1, 0, 0]
            #       [1, 0, 1, 0]
            #       [1, 0, 0, 1]
            
        elif mask_type == 'bidirectional_mlm':
            for i, tok in enumerate(token_to_decode):
                perm_mask[i, :, torch.arange(seq_len) != tok] = 0.0 # Allow full bidirectional context (masked-language-model-style. this is not autoregressive)
                # In query stream attention, tokens are always masked from themselves
                perm_mask[i, torch.arange(seq_len), torch.arange(seq_len)] = 1.0

                # [1, 0, 1]
                # [0, 1, 1]
                # [0, 0, 1]

                # [1, 1, 0]
                # [0, 1, 0]
                # [0, 1, 1]
                
                # [1, 0, 0]
                # [1, 1, 0]
                # [1, 0, 1]

        target_mapping = torch.zeros(
            (n_tokens, 1, seq_len), dtype=torch.float
        )  # Shape [batch_size=n_tokens, num_tokens_to_predict=1, seq_length=n_tokens]
        for i, tok in enumerate(token_to_decode):
            target_mapping[i, 0, tok] = 1.0 # Predict token tok at batch position i
            
        with torch.inference_mode():
            out = self.model(input_ids.to(self.device), perm_mask=perm_mask.to(self.device), target_mapping=target_mapping.to(self.device))
            # logits has shape [batch_size, 1, config.vocab_size] (1 is num_tokens_to_predict)
            index_corrected_logits = out.logits[:, 0, self.canonical_idx_to_xlnet_idx]
            # Ignore the last entry, corresponding to 'X'
            index_corrected_logits = index_corrected_logits[:, :-1]
            # Get temperature-adjusted probabilities
            probs = torch.nn.functional.softmax(index_corrected_logits / temperature, dim=-1)
            probs = probs / probs.sum(dim=-1, keepdim=True)

        return probs
        
class ProteinMPNNWrapper(nn.Module):
    def __init__(self, device=None):
        super().__init__()
        if device is not None:
            self.device = device
        elif torch.cuda.is_available():
            self.device = torch.device("cuda:0")
        else:
            self.device = torch.device('cpu')

        backbone_noise=0.00               # Standard deviation of Gaussian noise to add to backbone atoms
        #v_48_030=version with 48 edges 0.30A noise
        checkpoint_path ='./bayes_design/protein_mpnn/vanilla_model_weights/v_48_030.pt'
        hidden_dim = 128
        num_layers = 3
        checkpoint = torch.load(checkpoint_path, map_location=self.device, weights_only=True) 
        print('Number of edges:', checkpoint['num_edges'])
        noise_level_print = checkpoint['noise_level']
        print(f'Training noise level: {noise_level_print}')
        self.model = ProteinMPNN(num_letters=21, node_features=hidden_dim, edge_features=hidden_dim, hidden_dim=hidden_dim, num_encoder_layers=num_layers, num_decoder_layers=num_layers, augment_eps=backbone_noise, k_neighbors=checkpoint['num_edges'])
        self.model.to(self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.eval()
        print("Model loaded")
    
    def forward(self, seq, struct, decode_order, token_to_decode, mask_type, temperature=1.0):
        """Accept an amino acid sequence and protein structure 
        coordinates, return class probabilities for the next token
        Args:
            seq (len N list of len L_seq str): a list of string representations of an amino 
                acid sequence.
            struct ((L x 4 x 3) torch.Tensor): batch_size x seq_length x 
                num_atoms x num_coordinates tensor
            decode_order (len L list): list of the order of indices to decode.
                This determines the values in the permutation mask. Each index 
                attends to all indices that occur previously in the decoding_order.
            token_to_decode (len N tensor): tensor of indices in the range [0, L-1] indicating which 
                token to decode next.
        Returns:
            probs ((N x 20) torch.Tensor): a vector of probabilities for the next
                token
        """

        N = len(token_to_decode)
        L = len(seq[0])
        assert L == struct.shape[0], "Sequence length must match the number of residues in the provided structure"

        # Convert amino acid character to index
        seq = [re.sub(r"-", "X", s) for s in seq]
        seq = torch.tensor([[AMINO_ACID_ORDER.index(aa) for aa in s] for s in seq]).to(self.device)
        with torch.no_grad():
            if mask_type != 'bidirectional_mlm':
                decode_order = torch.tensor(decode_order).expand(N, L)
            elif mask_type == 'bidirectional_mlm':
                decode_order = torch.tensor(np.array([np.append(np.delete(decode_order, decode_order.index(tok)).tolist(), tok) for tok in token_to_decode]))
            
            struct = struct.expand(N, *struct.shape).to(self.device)
            # Default values
            # `mask` masks positions that are missing structural information
            mask = torch.isfinite(torch.sum(struct,(2,3))).to(torch.float32)
            isnan = torch.isnan(struct)
            struct[isnan] = 0.
            chain_M = torch.ones(N, L).float().to(self.device)
            chain_encoding_all = torch.ones(N, L).float().to(self.device)
            residue_idx = torch.arange(L).expand(N, L).to(self.device)

            _, logits = self.model(X=struct, S=seq, mask=mask, chain_M=chain_M, residue_idx=residue_idx, chain_encoding_all=chain_encoding_all, use_input_decoding_order=True, randn=None, decoding_order=decode_order.to(self.device))
            # Get temperature-adjusted probabilities
            probs = torch.nn.functional.softmax(logits / temperature, dim=-1)
            # N x L x 20
            # Ignore last entry, corresponding to 'X'
            probs = probs[:, :, :-1]
            probs = probs / probs.sum(dim=-1, keepdim=True)
            
        # Note that unlike XLNet, ProteinMPNN gives probabilities for all residues, not just the one to decode. 
        # So here we extract the probabilities for the residue to decode.
        return probs[range(len(token_to_decode)), token_to_decode]
        # N x 20

class BayesDesign(nn.Module):
    def __init__(self, device=None, bayes_balance_factor=0.002, **kwargs):
        super().__init__()
        if device is not None:
            self.device = device
        elif torch.cuda.is_available():
            self.device = torch.device("cuda:0")
        else:
            self.device = torch.device('cpu')

        self.seq_model = XLNetWrapper(device=device, **kwargs)
        self.seq_struct_model = ProteinMPNNWrapper(device=device, **kwargs)

        self.bayes_balance_factor = bayes_balance_factor

    def forward(self, seq, struct, decode_order, token_to_decode, mask_type='bidirectional_autoregressive', temperature=1.0):
        p_seq = self.seq_model(seq=seq, decode_order=decode_order, token_to_decode=token_to_decode, mask_type=mask_type, temperature=temperature).clone()
        p_seq_struct = self.seq_struct_model(seq=seq, struct=struct, decode_order=decode_order, token_to_decode=token_to_decode, mask_type=mask_type, temperature=temperature).clone()

        # Add a "balance factor" so that we don't end up with large probability ratios at the tails of the distributions
        p_seq += self.bayes_balance_factor
        p_seq_struct += self.bayes_balance_factor
        balanced_logits = (p_seq_struct / p_seq)
        
        # Normalize probabilities
        p_struct_seq = balanced_logits / balanced_logits.sum(dim=-1).unsqueeze(-1)

        return p_struct_seq


class TrRosettaWrapper():
    def __init__(self, data_location='./data/msa', database='./data/uniref30/UniRef30_2022_02', device=None):

        if device is not None:
            self.device = device
        elif torch.cuda.is_available():
            self.device = torch.device("cuda:0")
        else:
            self.device = torch.device('cpu')
        
        model_files = [*Path(DEFAULT_MODEL_PATH).glob('*.pt')]
        self.models = []
        for model_file in model_files[:1]:
            trRosetta = trRosettaNetwork(
                filters = 64,
                kernel = 3,
                num_layers = 61
            )
            trRosetta = trRosetta.to(self.device)
            trRosetta.load_state_dict(torch.load(model_file, map_location=self.device, weights_only=True))
            trRosetta.eval()
            self.models.append(trRosetta)

        self.data_location = data_location
        os.makedirs(self.data_location, exist_ok=True)

        self.database = database

    def __call__(self, seq, seq_id=""):
        """Pass a sequence through the trRosetta model and return the distogram
        Args:
            seq (str): a space-separated string representing the amino acid sequence
        Returns:
            distance ((L x L x 37) torch.Tensor): a distogram representing distance bin probabilities
        """
        seq_path = os.path.join(self.data_location, f'{seq_id}_{seq}.txt')
        seq_msa_path = os.path.join(self.data_location, f'msa_{seq_id}_{seq}.txt')
        # Make a fasta file with the sequence
        with open(seq_path, 'w') as f:
            f.write('>\n' + ''.join(seq.split()))
        # Get an MSA for the sequence
        if not os.path.exists(seq_msa_path):
            out = subprocess.run(['/root/hh-suite/bin/hhblits', '-v', '0', '-i', f'{seq_path}', '-oa3m', f'{seq_msa_path}', '-d', f'{self.database}'])
        x = preprocess(seq_msa_path).to(self.device)
        outputs = []
        with torch.no_grad():
            for model in self.models:
                output = model(x)
                outputs.append(output)
            averaged_outputs = [torch.stack(model_output).mean(dim=0).cpu().numpy().squeeze(0).transpose(1,2,0) for model_output in zip(*outputs)]
            # prob_theta, prob_phi, prob_distance, prob_omega
            output_dict = dict(zip(['theta', 'phi', 'dist', 'omega'], averaged_outputs))
            distance = output_dict['dist']
            distance = distance.squeeze()

        return distance

def test_xlnet_wrapper():
    import random

    seq, struct = get_protein()
    seq = seq[:10]
    given_characters_up_to = 5
    seq = ''.join([char if i < given_characters_up_to else '-' for i, char in enumerate(seq)])

    xlnet = XLNetWrapper()

    decoding_order = np.arange(len(seq)).tolist()
    random.shuffle(decoding_order)
    out_probs = xlnet(seq, decoding_order=decoding_order, token_to_decode=given_characters_up_to)

def test_protein_mpnn_wrapper():
    seq, struct = get_protein()
    given_characters_up_to = 20
    seq = ''.join([char if i < given_characters_up_to else '-' for i, char in enumerate(seq)])

    mpnn = ProteinMPNNWrapper()
    decoding_order = np.arange(len(seq)).tolist()
    out_probs = mpnn(seq=seq, struct=struct, decoding_order=decoding_order, token_to_decode=given_characters_up_to)



class TransformerXLU100Wrapper(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self):
        pass

class ProGenWrapper(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self):
        pass

class ESMIF1Wrapper(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self):
        pass

class PSSM():
    def __init__(self, pssm_path):
        # load pssm from pickle
        with open(pssm_path, 'rb') as f:
            self.pssm = pkl.load(f)

    def __call__(self, seq, struct, decode_order, token_to_decode, mask_type):
        return self.pssm[token_to_decode, :]


model_dict = {'xlnet':XLNetWrapper, 'protein_mpnn':ProteinMPNNWrapper, 'bayes_design':BayesDesign, 'pssm':PSSM, 'trRosetta':TrRosettaWrapper}
