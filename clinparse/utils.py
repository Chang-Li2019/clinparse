from varcode import Variant
from pyensembl import ensembl_grch38
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import gzip, io, os, re

def read_vcf(path):
    with gzip.open(path, 'r') as f:
        lines = [l.decode() for l in f if not l.decode().startswith('##')]
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
            ).rename(columns={'#CHROM': 'CHROM'})

def parse_clinvar(INFO):
    INFO = INFO.split(';')
    info_dict = {x.split('=')[0]:x.split('=')[1] for x in INFO}
    return info_dict.get('CLNSIG',None)

def get_top_effect(contig, start, ref, alt):
    myVariant = Variant(contig=contig, start=start, ref=ref, alt=alt, ensembl=ensembl_grch38)
    me = myVariant.effects()
    return me.top_priority_effect()

def is_protein_changing(chrom, pos, ref, alt):
    """
    Translate ClinVar genomic variants to proteomic variants
    """
    e = get_top_effect(chrom, pos, ref, alt)
    conseq = type(e).__name__
    if conseq == 'Substitution':
        return [e.gene_name, e.aa_mutation_start_offset, e.aa_alt]
    else:
        return [np.nan,np.nan,np.nan]

def process_clinvar(gz_file, use_high_quality_only=False):

    """
    Process ClinVar file in vcf.gz format.
    --------------------------------------
    Parameters:
        gz_file: input path of clinvar file
        use_high_quality_only: if set True, only "Pathogenic" and "Benign" variants will be included; 
                               otherwise, likely pathogenic/benign ones will also be included.
    """

    PATHOGENIC_LIST = ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]
    BENIGN_LIST = ["Benign", "Likely_benign", "Benign/Likely_benign"]
    if use_high_quality_only:
        PATHOGENIC_LIST = [PATHOGENIC_LIST[0]]
        BENIGN_LIST = [BENIGN_LIST[0]]

    clinvar_df = read_vcf(gz_file) # Read ClinVar as pandas dataframe
    ref_idx = clinvar_df['REF'].apply(lambda x: True if len(x)==1 else False) # Get index for single nucleotide variants
    alt_idx = clinvar_df['ALT'].apply(lambda x: True if len(x)==1 else False) # Get index for single nucleotide variants
    nucleotide_idx = (clinvar_df['REF']!='N') & (clinvar_df['ALT']!='N') # Get index for known nucleotide (remove 'N')
    pathogenic_conseq= clinvar_df['INFO'].apply(parse_clinvar)  # Get pathogenicity of variants
    pathogenic_idx = pathogenic_conseq.isin(PATHOGENIC_LIST) # get variants for one category
    benign_idx = pathogenic_conseq.isin(BENIGN_LIST)
    clinvar_tp = clinvar_df[ref_idx&alt_idx&nucleotide_idx&pathogenic_idx]  # Apply all previous filters to get qualifying variants
    clinvar_tn = clinvar_df[ref_idx&alt_idx&nucleotide_idx&benign_idx] 
    ns_data_tp = clinvar_tp.iloc[:,[0,1,3,4]].apply(lambda x: is_protein_changing(x['CHROM'], x['POS'], x['REF'], x['ALT']), axis=1) # Get protein changing consequences.
    ns_data_tn = clinvar_tn.iloc[:,[0,1,3,4]].apply(lambda x: is_protein_changing(x['CHROM'], x['POS'], x['REF'], x['ALT']), axis=1) # Get protein changing consequences.
    clinvar_ns_tp = pd.DataFrame(ns_data_tp.to_list()) # process data to create new data frame
    clinvar_ns_tp = clinvar_ns_tp[-clinvar_ns_tp[0].isna()] # remove missing values

    clinvar_ns_tn = pd.DataFrame(ns_data_tn.to_list()) # process data to create new data frame
    clinvar_ns_tn = clinvar_ns_tn[-clinvar_ns_tn[0].isna()] # remove missing values

    clinvar_ns_tp.columns = ['gene_name', 'aa_pos', 'aa_alt']  # rename columns
    clinvar_ns_tn.columns = ['gene_name', 'aa_pos', 'aa_alt']
    clinvar_ns_tp['label'] = 1
    clinvar_ns_tn['label'] = 0
    clinvar_ns_tp['aa_pos'] = clinvar_ns_tp['aa_pos'].astype(np.int32)+1  # change the types of the column and compensate for 0-based index
    clinvar_ns_tn['aa_pos'] = clinvar_ns_tn['aa_pos'].astype(np.int32)+1
    clinvar_ns = pd.concat([clinvar_ns_tp, clinvar_ns_tn])
    return clinvar_ns


def parse_protein_seq(file_path = './data/Ensembl107ProteinSeq.fasta'):
    from Bio import SeqIO
    x = list(SeqIO.parse(file_path, "fasta"))
    return {record.id.split('|')[1]:record.seq._data.decode().replace('*','') for record in x}


def read_pdb(file):
     with open(file,'r') as f:
          d = f.readlines()
     d = [x.split() for x in d]
     d = [x for x in d if ((len(x)>2) and (x[2]=='CA'))]
     plddt = [float(x[-2]) for x in d if x[0]=='ATOM']
     coords = [x[6:9] for x in d if x[0]=='ATOM']
     correct_coords = []
     for coord in coords:  ## TODO:more efficient way to do this?
          correct_coord = []
          for x in coord:
               correct_coord.extend(correct_pdb_line(x))
          correct_coord = correct_coord if len(correct_coord)==3 else correct_coord[:-2] if len(correct_coord)==5 else correct_coord[:-1]
          correct_coords.append(correct_coord)
          
     df1 = pd.DataFrame(correct_coords, columns=['x','y','z'])
     df2 = pd.DataFrame(plddt, columns=['plDDT'])
     return pd.concat([df1,df2], axis=1).astype('float')
     
def correct_pdb_line(coord):
     """
     Some pdb files from AlphaFold have two or three un-delimited coordinates
     """
     match = list(re.finditer('-[1-9]', coord))
     if match==[]:
          return [float(coord)]
     else: 
          x = [i.start() for i in match if i.start()!=0]
          if x==[]:
               return [float(coord)]
          elif len(x)==1:
               return [float(coord[:x[0]]), float(coord[x[0]:])]
          else:
               return [float(coord[:x[0]]), float(coord[x[0]:x[1]]), float(coord[x[1]:])]


def rmsd(V, W):
     """
     rmsd between two 1/2d arrays
     """
     diff = np.array(V) - np.array(W)
     N = len(V)
     return np.sqrt((diff * diff).sum() / N)

def plot_contact(pdf_df,name='ref'):
    mat = cdist(pdf_df[['x','y','z']],pdf_df[['x','y','z']])
    # bins = np.histogram(mat, bins = 6)[1]  # as there was no global optimal cutoff exists (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1983-3)
    bins = np.array([0,6,12,20,10000])
    digitized = np.digitize(mat, bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(digitized, cmap='Greys_r')
    fig.colorbar(cax)
    fig.savefig(f'{name}.pdf')