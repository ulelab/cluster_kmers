import pandas as pd
from skbio import RNA
from skbio.alignment import global_pairwise_align_nucleotide
from itertools import combinations
import warnings
import logomaker
import matplotlib.pyplot as plt
import seaborn as sns
import seqlogo
warnings.filterwarnings("ignore", message="You're using skbio's python implementation of Needleman-Wunsch alignment.")

# To be able to access IUPAC codes from this script
# Get the directory of the current file (script C)
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
# Construct the path to the resource file
resource_file_path = os.path.join(script_dir, '../res/', 'IUPAC_nt_code.csv')


def get_prefix_suffix(j, motif):
    prefix = j.split(motif)[0]
    suffix = j.split(motif)[-1]
    return prefix, suffix

def get_alignments(motif_pairs_list, match_s, missm_s):
    df_cons = pd.DataFrame(columns=['consensus', 'score', 'aligned_mots', 'motif_pair'])
    for i, pair in enumerate(motif_pairs_list):
        alignment, score, start_end_positions = global_pairwise_align_nucleotide(
            RNA(pair[0]), RNA(pair[1]),
            match_score=match_s, mismatch_score=missm_s
            )
        df_cons.loc[i, :] = str(alignment.consensus()), score, tuple(str(alignment).split('\n')[-2:]), pair
        df_cons = df_cons.sort_values(by='score', ascending=False)
    return df_cons

def get_alignment(mots, match_s=2, missm_s=-1):
    motif_pairs = list(combinations(mots, r=2))
    motif_pairs = [tuple(sorted(m)) for m in motif_pairs]
    #Get pairwise alignments
    df_cons = get_alignments(motif_pairs, match_s, missm_s)
    # Get best alignment out of all as a tuple of aligned motifs
    aligned = list(df_cons.sort_values(by='score', ascending=False)['motif_pair'].values.tolist()[0])
    aligned_sequences = list(df_cons.sort_values(by='score', ascending=False)['aligned_mots'].values.tolist()[0])
    # remove aligned motifs from the list of all motifs
    remaining_mots = [m for m in mots if m not in aligned]
    while len(remaining_mots)>0:
        #Find the best alignment related to one of the aligned motifs (reference motifs)
        aligned_pairs = list(combinations(aligned, r=2))
        aligned_pairs = [tuple(sorted(m)) for m in aligned_pairs]
        lookup_pairs = [m for m in motif_pairs if m not in aligned_pairs]
        lookup_pairs = [m for m in lookup_pairs if (m[0] in aligned) or (m[1] in aligned)]
        df_select = df_cons.loc[df_cons['motif_pair'].isin(lookup_pairs)].sort_values(by='score', ascending=False)
        second_alignment = df_select['aligned_mots'].values.tolist()[0]
        second_pair = df_select['motif_pair'].values.tolist()[0]
        # Set reference motif
        y = [m for m in aligned if m in second_pair][0]
        ref1 = [m for m in aligned_sequences if m.replace('-','')==y][0]
        ref2 = [m for m in second_alignment if m.replace('-','')==y][0]
        # Define other motif to align
        seq2 = [m for m in second_alignment if m.replace('-','')!=y][0]
        # All other aligned motifs must be changed accordingly
        al_seqs = [m for m in aligned_sequences if m.replace('-', '')!=y]
        #The ref2 and seq2 will get this prefix and suffix
        p2, s2 = get_prefix_suffix(ref1, y)
        #The ref1 and seq1 will get this prefix and suffix
        p1, s1 = get_prefix_suffix(ref2, y)
        #New sequences
        ref1 = p1 + ref1 + s1
        al_seqs = [p1 + seq1 + s1 for seq1 in al_seqs]
        seq2 = p2 + seq2 + s2
        aligned_sequences = al_seqs + [ref1, seq2]
        aligned = [s.replace('-', '') for s in aligned_sequences]
        remaining_mots = [m for m in mots if m not in aligned]
    return aligned_sequences

def get_consensus(aligned, df_iupac_codes=resource_file_path, l=5):
    #Generate consensus
    df_consensus = pd.DataFrame(0, columns=[i for i in range(len(aligned[0]))], index=['U', 'G', 'C', 'A'])
    for s in aligned:
        for pos, letter in enumerate(s):
            if letter in df_consensus.index.tolist():
                df_consensus.loc[letter, pos] += 1
    #Make a full scored consensus
    df_sequence = pd.DataFrame(index=df_consensus.columns, columns=['nt', 'score'])
    for pos in df_consensus.columns.tolist():
        s = df_consensus[pos].where(df_consensus[pos] == df_consensus[pos].max()).dropna()
        #Only one max nucleotide
        if len(s)==1:
            df_sequence.loc[pos, 'nt'] = s.index.tolist()[0]
            df_sequence.loc[pos, 'score'] = s.values.tolist()[0]
        #If there is a tie - write both as sequence
        else:
            bases_tied = ', '.join(sorted(s.index.tolist()))
            code = df_iupac_codes.loc[df_iupac_codes['bases_sorted']==bases_tied, 'IUPAC_code'].values.tolist()[0]
            df_sequence.loc[pos, 'score'] = s.sum()
            df_sequence.loc[pos, 'nt'] = code
    #use a sliding window of kmer length to determine consensus sequence
    df_sequence['rolling_sum'] = df_sequence['score'].rolling(l).sum()
    end = df_sequence['rolling_sum'].idxmax()
    start = end - (l-1)
    consensus = ''.join(df_sequence.loc[start:end, 'nt'].values.tolist())
    return consensus, df_sequence

def get_consensus_from_motifs(mots, l=4, iupac_path=resource_file_path, match_s=2, missm_s=-1):
    aligned_sequences = get_alignment(mots, match_s, missm_s)
    #Import iupac code to make consensuses
    df_iupac_codes = pd.read_csv(iupac_path, sep='\t')
    df_iupac_codes['bases_sorted'] = [', '.join(sorted(s.split(' or '))) for s in df_iupac_codes['Base'].values.tolist()]
    consensus = get_consensus(aligned_sequences, df_iupac_codes, l)[0]
    return consensus, aligned_sequences

def GetPfm(aligned):
    """
    Input: aligned sequences from get_alignment.
    """
    # Generate position frequency matrix from alignment
    idx = [i for i in range(len(aligned[0]))]
    df_pfm = pd.DataFrame(0, index=idx, columns=['A', 'C', 'G', 'U', 'N', '-'])
    # print('Alignment:')
    for s in aligned:
        # print(s)
        for pos, letter in enumerate(s):
            df_pfm.loc[pos, letter] += 1
    # Remove positions in PFM, where no nucleotides
    df_pfm = df_pfm.loc[(df_pfm[['A', 'C', 'G', 'U']] != 0).any(axis=1)]
    df_pfm.reset_index(drop=True, inplace=True)
    return df_pfm

def PlotKmerLogoLogomaker(kmer_list, match_s=2, missm_s=-1, colorScheme='default_RNA'):
    if colorScheme == 'default_RNA':
        colorScheme = {
        'A': '#48AD7E',
        'U': '#FC5252',
        'G': '#F5D741',
        'C': '#55BDF5',
        '-': 'white',
        'N': 'white',
        }
    if len(kmer_list) > 1:
        aligned = get_alignment(kmer_list, match_s, missm_s)
    elif len(kmer_list) == 1:
        aligned = kmer_list
    df_pfm = GetPfm(aligned)
    # convert pfm to ppm and pwm
    pfm = seqlogo.Pfm(df_pfm, alphabet_type='reduced RNA', background=0)
    # convert pfm to position probability matrix
    # ppm = cpm.ppm
    ppm = seqlogo.pfm2ppm(pfm, alphabet_type='reduced RNA')
    nts = ['A', 'U', 'G', 'C']
    ppm = ppm[nts]
    fig, ax = plt.subplots()
    logomaker.Logo(ppm, color_scheme=colorScheme, ax=ax)
    return fig, ax

