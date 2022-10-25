#!/usr/bin/env python3
"""
Author : changli <changli@localhost>
Date   : 2022-10-24
Purpose: Given an input gene name, get missense pathogenic/benign variants in AminoAcid positions
"""

import argparse, os
from typing import NamedTuple, TextIO
import pandas as pd

from utils import process_clinvar, parse_protein_seq

class Args(NamedTuple):
    """ Command-line arguments """
    positional: str
    file: TextIO
    outfile: TextIO

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='clinvar parser',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('positional',
                        metavar='str',
                        help='Gene name of interest',
                        default=None)

    parser.add_argument('-f',
                        '--file',
                        help='ClinVar vcf file',
                        metavar='FILE',
                        type=argparse.FileType('rb'),
                        default=None)

    parser.add_argument('-o',
                        '--outfile',
                        help='Output file name',
                        metavar='FILE',
                        type=argparse.FileType('w'),
                        default="cand_gene_output")


    args = parser.parse_args()

    return Args(args.positional, args.file, args.outfile)


# --------------------------------------------------
def main() -> None:

    args = get_args()
    file_arg = args.file
    gene_name = args.positional
    outfile_arg = args.outfile

    print('ClinVar source file = "{}"'.format(file_arg.name if file_arg else ''))
    print('Output file = "{}"'.format(outfile_arg.name if outfile_arg else ''))
    print(f'Gene name of interest = "{gene_name}"')
    if os.path.exists('ClinVar.db'):
        clinvar_ns = pd.read_csv('ClinVar.db', sep=',')
    else:
        print('Processing ClinVar variants. This only needs to run once, which may take around 20-30 minutes ...')
        clinvar_ns = process_clinvar(file_arg)
        clinvar_ns.to_csv('ClinVar.db', sep=',')

    if gene_name is not None:
        idx = (clinvar_ns['gene_name'] == gene_name)
        if idx.any():
            clinvar_ns[idx].to_csv(outfile_arg, sep='\t', index=False)
            gene_to_seq = parse_protein_seq()
            gene_to_seq[gene_name]
            if not os.path.exists(f'{gene_name}_mutant'):
                os.makedirs(f'{gene_name}_mutant')
            refseq = gene_to_seq.get(gene_name, None)
            if refseq is not None:
                for _,rowData in clinvar_ns[idx].iterrows():
                    alt = rowData['aa_alt']
                    pos = rowData['aa_pos']
                    label = rowData['label']
                    filename = f"{gene_name}_{pos}_{alt}_{label}"
                    outstring = f">{gene_name}_{pos}_{alt}_{label}\n"
                    altseq = list(refseq)
                    altseq[pos-1] = alt
                    altseq = ''.join(altseq)+'\n'
                    outstring += altseq
                    with open(os.path.join(f"{gene_name}_mutant", filename), 'w') as f:
                        f.write(outstring)
        else:
            print("No gene found for name: {gene_name}")
        
       
    
# --------------------------------------------------
if __name__ == '__main__':
    main()
