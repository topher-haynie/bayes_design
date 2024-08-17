# scRMSD vs. ProteinMPNN in silico experiment

from Bio.PDB import PDBList, PDBParser
import torch
import numpy as np
import pandas as pd
import os
import sys
import argparse
from bayes_design.decode import decode_order_dict, decode_algorithm_dict
from bayes_design.model import model_dict
from bayes_design.utils import get_protein, align_and_crop, get_ball_mask, get_fixed_position_mask

from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import PPBuilder


def get_sequence(chain):
    """Extract the sequence of a chain."""
    ppb = PPBuilder()
    for pp in ppb.build_peptides(chain):
        return pp.get_sequence()

def align_sequences(seq1, seq2):
    """Perform global alignment between two sequences."""
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = max(alignments, key=lambda x: x[2])  # Get alignment with highest score
    return best_alignment

def calculate_identity(align1, align2):
    """Calculate sequence identity between two aligned sequences."""
    matches = sum(res1 == res2 for res1, res2 in zip(align1, align2))
    return matches / len(align1) * 100

def find_pdb_chain_matches():
    """For each protein in the PDB calculate the scRMSD between any chain in that protein and any chains 

    Pseudocode:

    # Iterate over the PDB
    # For each protein
        # For each chain in that protein
            # Iterate over the PDB
            matches = dict()
            # For each chain in each protein
                # If there is chain has 98% sequence id to that chain*
                    # Create variables identifying the beginning residue id and end residue id of the overlapping portions between the chains
                    # Calculate the residue-wise RMSD between the two chains under a kabsch alignment
                    # Iterate over each 10-residue linear span and calculate the average RMSD for that span.
                    # Identify the top linear RMSD span
                    # Create a dict with this info:
                    matches[new_pdb_id + "_" + new_chain_id]["this_pdb_residue_ids"] = (this_pdb_beginning_residue_id, this_pdb_end_residue_id)
                    matches[new_pdb_id + "_" + new_chain_id]["new_pdb_residue_ids"] = (beginning_residue_id, end_residue_id)
                    matches[new_pdb_id + "_" + new_chain_id]["rmsds"] = rmsd
                    matches[new_pdb_id + "_" + new_chain_id]["linear_rmsds"] = linear_rmsds
                    matches[new_pdb_id + "_" + new_chain_id]["top_linear_rmsd"] = top_linear_rmsd
    # *98% sequence identity means that one chain is a 98% sequence identity match to the other chain. It may cover 60% of the residues in the chain, but the 60% that are covered are 98% identical to the other chain. 
    """
    pdbl = PDBList()
    all_pdb_ids = pdbl.get_all_entries()

    for pdb_id in all_pdb_ids:
        pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
        parser = PDBParser()
        structure = parser.get_structure(pdb_id, f"{pdb_id}.pdb")
        for model in structure:
            for chain in model:
                seq1 = get_sequence(chain)
                matches = dict()
                for new_pdb_id in all_pdb_ids:
                    if new_pdb_id == pdb_id:
                        continue  # Skip comparing the same structure
                    pdbl.retrieve_pdb_file(new_pdb_id, pdir='.', file_format='pdb')
                    new_structure = parser.get_structure(new_pdb_id, f"{new_pdb_id}.pdb")
                    for new_model in new_structure:
                        for new_chain in new_model:
                            seq2 = get_sequence(new_chain)
                            alignment = align_sequences(seq1, seq2)
                            identity = calculate_identity(alignment[0], alignment[1])

                            if identity >= 95:
                                start_res1 = alignment[3]
                                end_res1 = alignment[4]
                                start_res2 = alignment[5]
                                end_res2 = alignment[6]
                                matches[(pdb_id, chain.id)] = {
                                    "matching_pdb_id": new_pdb_id,
                                    "identity": identity,
                                    "overlap_residue_range_1": (start_res1, end_res1),
                                    "overlap_residue_range_2": (start_res2, end_res2),
                                }
                                print(f"Match found: {pdb_id}:{chain.id} with {new_pdb_id}:{new_chain.id}")
                                print(f"Identity: {identity:.2f}%")
                                print(f"Overlap Residues: Chain1 ({start_res1}, {end_res1}), Chain2 ({start_res2}, {end_res2})")
             
                                






"""Iterate over matches and identify the PDB chains corresponding to the top 100 'top_linear_rmsd' values. Inverse fold these chains with ProteinMPNN and with CSDesign. Remove duplicates (if a matches b, then there will be an entry for b matches a)."""

