#!/usr/bin/env python3
"""
Author : changli <changli@localhost>
Date   : 2022-11-01
Purpose: Rock the Casbah
"""

import argparse, os, re
from typing import NamedTuple, TextIO
import matplotlib.pyplot as plt


import pandas as pd
import numpy as np

from utils import read_pdb, plot_contact, rmsd

class Args(NamedTuple):
    """ Command-line arguments """
    tpfolder: str
    tnfolder: str
    file: str

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def file_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)
# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """
    parser = argparse.ArgumentParser(
        description='Rock the Casbah',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-tp',
                        '--tpfolder',
                        help='.pdf files for pathogenic variants',
                        metavar='str',
                        type=dir_path,
                        default=None)

    parser.add_argument('-tn',
                        '--tnfolder',
                        help='.pdf files for benignvariants',
                        metavar='str',
                        type=dir_path,
                        default=None)

    parser.add_argument('-f',
                        '--file',
                        help='A reference pdb file for the gene of interest',
                        metavar='FILE',
                        type=file_path)


    args = parser.parse_args()

    return Args(args.tpfolder, args.tnfolder, args.file)


# --------------------------------------------------
def main() -> None:
    """ Visualize the plDDT vs position """

    args = get_args()
    tpfolder = args.tpfolder
    tnfolder = args.tnfolder
    ref_file = args.file

    print(f'TP folder path = {tpfolder}')
    print(f'TN folder path = {tnfolder}')
    print(f'referece file path = {ref_file}')

    ref_df = read_pdb(ref_file)
    tp_df = {}
    tn_df = {}
    if tpfolder is not None:
        for file in os.listdir(tpfolder):
            ### only taking relaxed model right now
            if file.endswith('_relaxed_rank_1_model_3.pdb'):
                df = read_pdb(os.path.join(tpfolder, file))
                idx = list(re.finditer('_',file))
                gene_idx = idx[0].start()
                gene = file[:gene_idx]
                key = file[idx[0].end():idx[2].start()]
                if gene not in tp_df:
                    tp_df[gene] = {key:df}
                else:
                    tp_df[gene].update({key:df})
        print(f'{len(tp_df[gene])} positive variants processed.')
    if tnfolder is not None:
        for file in os.listdir(tnfolder):
            if file.endswith('_relaxed_rank_1_model_3.pdb'):
                df = read_pdb(os.path.join(tnfolder, file))
                idx = list(re.finditer('_',file))
                gene_idx = idx[0].start()
                gene = file[:gene_idx]
                key = file[idx[0].end():idx[2].start()]
                if gene not in tn_df:
                    tn_df[gene] = {key:df}
                else:
                    tn_df[gene].update({key:df})
        print(f'{len(tn_df[gene])} negative variants processed.')

    #### plddt vs amino acid position
    #### Contact map


    fig, ax = plt.subplots(1)
    ax.scatter(ref_df.index+1, ref_df.plDDT, color='grey', alpha=0.2, label='All variants')
    
    if tp_df != {}:
        tp_rmsd = {}
        tp_loci = []
        for k,v in tp_df[gene].items():
            pos = k.split('_')[0]
            plot_contact(v, name=f'{k}')
            tp_loci.append(int(pos)-1)
            tp_rmsd[k] = rmsd(v[['x','y','z']], ref_df[['x','y','z']])
        tp_plddt = ref_df[ref_df.index.isin(tp_loci)]
        ax.scatter(tp_plddt.index+1, tp_plddt.plDDT, color='r', label='Pathogenic')
    if tn_df != {}:
        tn_rmsd = {}
        tn_loci = []
        for k,v in tn_df[gene].items():
            pos = k.split('_')[0]
            plot_contact(v, name=f'{k}')
            tn_loci.append(int(pos)-1)
            tn_rmsd[k] = rmsd(v[['x','y','z']], ref_df[['x','y','z']])
        tn_plddt = ref_df[ref_df.index.isin(tn_loci)]
        ax.scatter(tn_plddt.index+1, tn_plddt.plDDT, color='g', label='Benign')
    ax.set_xlabel('Amino acid sequence')
    ax.set_ylabel('plDDT values')
    fig.legend()
    fig.tight_layout()
    fig.savefig(f'{gene}_plDDT.pdf')
    plot_contact(ref_df)

    ### save rmsd values
    if tn_rmsd != {}:
        pd.DataFrame.from_dict(tn_rmsd, orient='index').to_csv(f'tn_rmsd_{gene}.tsv', sep='\t')
    if tp_rmsd != {}:
        pd.DataFrame.from_dict(tp_rmsd, orient='index').to_csv(f'tp_rmsd_{gene}.tsv', sep='\t')

# --------------------------------------------------
if __name__ == '__main__':
    main()
