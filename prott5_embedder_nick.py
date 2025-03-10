#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:33:22 2020

@author: mheinzinger

Refactored 2025-03-07 by Nick Pullen @ univie
"""

import argparse
import time
import os
import logging
from pathlib import Path
import sys

import torch
import h5py
from transformers import T5EncoderModel, T5Tokenizer

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
print("Using device: {}".format(device))

def get_T5_model(model_dir, transformer_link = "Rostlab/prot_t5_xl_half_uniref50-enc"):
    print("Loading: {}".format(transformer_link))
    if model_dir is not None:
        print("##########################")
        print("Loading cached model from: {}".format(model_dir))
        print("##########################")
    model = T5EncoderModel.from_pretrained(transformer_link, cache_dir=model_dir)
    # only cast to full-precision if no GPU is available
    if device==torch.device("cpu"):
        print("Casting model to full precision for running on CPU ...")
        model.to(torch.float32)

    model = model.to(device)
    model = model.eval()
    vocab = T5Tokenizer.from_pretrained(transformer_link, do_lower_case=False )
    return model, vocab


def read_fasta( fasta_path ):
    '''
        Reads in fasta file containing multiple sequences.
        Returns dictionary of holding multiple sequences or only single 
        sequence, depending on input file.
    '''
    
    sequences = dict()
    with open( fasta_path, 'r' ) as fasta_f:
        for line in fasta_f:
            # get uniprot ID from header and create new entry
            if line.startswith('>'):
                uniprot_id = line.replace('>', '').strip()
                # replace tokens that are mis-interpreted when loading h5
                uniprot_id = uniprot_id.replace("/","_").replace(".","_")
                sequences[ uniprot_id ] = ''
            else:
                # repl. all white-space chars and join seqs spanning multiple lines
                sequences[ uniprot_id ] += ''.join( line.split() ).upper().replace("-","") # drop gaps and cast to upper-case
                
    return sequences


def get_embeddings(seq_path, 
                   emb_path, 
                   model_dir, 
                   per_protein, # whether to derive per-protein (mean-pooled) embeddings
                   max_residues=4000, # number of cumulative residues per batch
                   max_seq_len=1000, # max length after which we switch to single-sequence processing to avoid OOM
                   max_batch=100 # max number of sequences per single batch
                   ):
    
    seq_dict = dict()

    # Read in fasta
    seq_dict = read_fasta( seq_path )
    model, vocab = get_T5_model(model_dir)

#    print('########################################')
#    print('Example sequence: {}\n{}'.format( next(iter(
#            seq_dict.keys())), next(iter(seq_dict.values()))) )
#    print('########################################')
    print('Total number of sequences: {}'.format(len(seq_dict)))

    avg_length = sum([ len(seq) for _, seq in seq_dict.items()]) / len(seq_dict)
    n_long     = sum([ 1 for _, seq in seq_dict.items() if len(seq)>max_seq_len])
    seq_dict   = sorted( seq_dict.items(), key=lambda kv: len( seq_dict[kv[0]] ), reverse=True )
    
    print("Average sequence length: {}".format(avg_length))
    print("Number of sequences >{}: {}".format(max_seq_len, n_long))

    # Here we flush all the prior print statements to stdout so we can check that sth is happening
    # Otherwise this gets held in a buffer until job completion which is not very helpful
    sys.stdout.flush()

    # Initialize counters
    batch_count = 0 # Batches for logging
    new_embeddings_count = 0  # How many new embeddings were processed

    start = time.time()
    batch = list()
    for seq_idx, (pdb_id, seq) in enumerate(seq_dict,1):
        seq = seq.replace('U','X').replace('Z','X').replace('O','X')
        seq_len = len(seq)
        seq = ' '.join(list(seq))
        batch.append((pdb_id,seq,seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed 
        n_res_batch = sum([ s_len for  _, _, s_len in batch ]) + seq_len 
        if len(batch) >= max_batch or n_res_batch>=max_residues or seq_idx==len(seq_dict) or seq_len>max_seq_len:
            # Checkpointing - Open the H5 file to get the already processed IDs
            try:
                with h5py.File(str(emb_path), "r") as hf:
                    processed_ids = set(hf.keys())
            except OSError:
                # File doesn't exist yet, so no IDs have been processed
                processed_ids = set()

            # Unpack the current batch
            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()
            batch_count += 1

            # Filter out sequences that have already been processed
            to_process = [(pid, seq, s_len) for pid, seq, s_len in zip(pdb_ids, seqs, seq_lens) if pid not in processed_ids]
    
            # Calculate total batch length using all sequences in the batch (for logging)
            # should be lower than max_residues
            total_batch_length = sum(seq_lens)
            all_ids_status = "\n".join(
                f"Batch {batch_count}: {'NEW' if pid not in processed_ids else 'EXISTING'} - {pid} (L={s_len})"
                for pid, s_len in zip(pdb_ids, seq_lens)
            )

            log_message = (
                f"Batch {batch_count}: Total batch length: {total_batch_length}, "
                f"{len(to_process)} new sequences, {len(pdb_ids) - len(to_process)} previous sequences.\n"
                f"IDs:\n{all_ids_status}\n"
            )

            # If no new sequences need processing, log the batch summary anyway and skip any computation
            if not to_process:
                logging.info(log_message)
                continue
    
            # These are the new sequences that need processing
            proc_ids, proc_seqs, proc_seq_lens = zip(*to_process)

            token_encoding = vocab.batch_encode_plus(proc_seqs, add_special_tokens=True, padding="longest")
            input_ids      = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)
            
            try:
                with torch.no_grad():
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
            except RuntimeError:
                # We record which batch failed and (all) its constituent proteins
                all_ids_status_fail = "\n".join(
                    f"Batch {batch_count}: FAIL - {pid} (L={s_len})"
                    for pid, s_len in zip(proc_ids, proc_seq_lens)
                )
                log_message = (
                    f"Batch {batch_count}: Total batch length: {total_batch_length}, "
                    f"{len(to_process)} new sequences, {len(pdb_ids) - len(to_process)} previous sequences.\n"
                    f"IDs:\n{all_ids_status_fail}\n"
                    f"FAIL: Batch {batch_count} encountered RuntimeError and will be skipped.\n"
                )
                # Save FAILs to log as well so that all proteins will have either NEW, EXISTING or FAIL in the log file
                logging.info(log_message)
                # This will go to the .out file and should indicate the exact protein that failed in the batch
                print("Batch {} with total batch length {} RuntimeError during embedding for {} (Length={} AAs). Try lowering batch size. ".format(batch_count, total_batch_length, proc_ids[-1], proc_seq_lens[-1]) +
                      "If single sequence processing does not work, you need more vRAM to process your protein.")
                sys.stdout.flush()
                continue

            # batch-size x seq_len x embedding_dim
            # extra token is added at the end of the seq
            with h5py.File(str(emb_path), "a") as hf:
                for batch_idx, identifier in enumerate(proc_ids):
                    s_len = proc_seq_lens[batch_idx]
                    # slice-off padded/special tokens
                    emb = embedding_repr.last_hidden_state[batch_idx,:s_len]
                    
                    if per_protein:
                        emb = emb.mean(dim=0)
                
                    hf.create_dataset(identifier, data=emb.detach().cpu().numpy().squeeze())
                    new_embeddings_count += 1
                    #print("new_embeddings_count=", new_embeddings_count)

            # Append the completion details to the log message after processing each batch
            log_message += f"Completed batch {batch_count}: Processed {len(proc_ids)} new sequences.\n"
            logging.info(log_message)

    end = time.time()

    print('\n############# OVERALL STATS #############')
    print('Total new embeddings processed in this run: {}'.format(new_embeddings_count))
    print('Total time: {:.2f}[s]; time/prot: {:.4f}[s]; avg. len of all proteins= {:.2f}'.format( 
            end-start, (end-start)/new_embeddings_count if new_embeddings_count > 0 else 0, avg_length))
    return True


def create_arg_parser():
    """Creates and returns the ArgumentParser object."""

    # Instantiate the parser
    parser = argparse.ArgumentParser(description=( 
            't5_embedder.py creates T5 embeddings for a given text '+
            ' file containing sequence(s) in FASTA-format.') )
    
    # Required positional argument
    parser.add_argument( '-i', '--input', required=True, type=str,
                    help='A path to a fasta-formatted text file containing protein sequence(s).')

    # Required positional argument
    parser.add_argument( '-o', '--output', required=True, type=str, 
                    help='A path for saving the created embeddings as NumPy npz file.')

    # Optional positional argument
    parser.add_argument('--model', required=False, type=str,
                    default=None,
                    help='A path to a directory holding the checkpoint for a pre-trained model' )

    # Optional argument
    parser.add_argument('--per_protein', type=int, 
                    default=1,
                    help="Whether to return per-residue embeddings (0) or the mean-pooled per-protein representation (1: default). \
                          Changed from original!")

    # Optional positional argument
    parser.add_argument( '-l', '--log', required=False, type=str, 
                    help='A path for saving a logging file. Otherwise written to stdout')
    return parser

def main():
    parser     = create_arg_parser()
    args       = parser.parse_args()
    
    seq_path   = Path( args.input )
    emb_path   = Path( args.output)
    model_dir  = Path( args.model ) if args.model is not None else None

    per_protein = False if int(args.per_protein)==0 else True
    log_path = Path(args.log) if args.log is not None else None

    # Very minor, but if the log file exists, we open it and write a newline to separate the appearance of jobs
    if log_path is not None and os.path.exists(log_path):
        with open(log_path, "a") as f:
            f.write("****************************************************************************************************************************\n\n")
    logging.basicConfig(filename=log_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - [Process: %(process)d - %(message)s')
    
    get_embeddings( seq_path, emb_path, model_dir, per_protein=per_protein )

if __name__ == '__main__':
    print("Starting...")
    torch.cuda.empty_cache()
    print('Number of threads = {}'.format(torch.get_num_threads()))
    main()
    logging.info("=========DONE============")
# Example to run:
# python prott5_embedder_nick.py --input Ecoli/494lines.fasta --output Ecoli/protein_embeddings.h5 --log /lisc/scratch/cube/pullen/testing.log
