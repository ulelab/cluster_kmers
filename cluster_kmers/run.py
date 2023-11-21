import sklearn
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import zscore
from . import get_clustering as clust
from . import get_alignment as align
from distutils.util import strtobool
import argparse
import os

def cli():
    parser = argparse.ArgumentParser(
        prog = 'cluster_kmers',
        description='Cluster k-mers based on sequence or on a combination of sequence and enrichment in CLIP data.'
        )
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-k',"--kmers", type=str, nargs='+', required=True, dest='kmers',
                        help='A list of k-mers encoded in RNA alphabet separated by spaces: AAGG GGAG GCCU.')
    required.add_argument('-o',"--output_folder", type=str, required=True, dest='results_folder',
                        help='A path to an existing output folder, for example "~/results"')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-co',"--cluster_on", type=str, choices=['seq', 'seq_enrichment'], dest='cluster_on', default='seq', required=False,
                        help="""
                        Inputs to clustering. Valid options are:
                        seq - Cluster only based on sequence similarity.
                        seq_enrichment - Cluster based on sequence similarity and based on enrichment of motifs in CLIP data
                                         (this option requires arguments passed to -peka)
                        """)
    optional.add_argument('-n',"--n_clusters", dest='n_clusters', type=(lambda x: 'auto' if str(x) == 'auto' else int(x)), default='auto', required=False,
                        help='Number of clusters to split k-mers into. Valid options are "auto" or integer. Default is "auto".')
    optional.add_argument('-peka',"--peka_files", nargs='+', dest='peka_files', default=None, required=False,
                        help='A list of peka output files with extensions *mer_distribution_{region_name}.tsv, separated by spaces.')
    optional.add_argument('-tl',"--min_token_length", type=int, default=1, dest='min_token_length', required=False,
                        help="""
                        Minimal length of a substrings used for clustering. For k-mers with lengths greater than 5, setting this value
                        to be greater than 1 can improve the results of clustering.
                        """)
    optional.add_argument('-wl',"--weblogos", dest='weblogos', default='True', type=(lambda x:bool(strtobool(x))), choices=[True, False], required=False,
                        help='Whether to plot weblogos for motif groups, True by default.')
    optional.add_argument('-cons',"--consensus_length", dest='consensus_length', default='auto', type=(lambda x: 'auto' if str(x) == 'auto' else int(x)), required=False,
                        help="""
                        Length of consensus sequence to name k-mer groups. Automatically this length is determined as k-mer length - 1.
                        Valid choices are "auto" or integer.
                        """)
    args = parser.parse_args()
    print(args)

    return args.kmers, args.cluster_on, args.n_clusters, args.peka_files, args.results_folder, args.min_token_length, args.weblogos, args.consensus_length


def GetPEKAScoresFromTsv(file_list):
    """
    Combines PEKA-scores columns from a list of tsv files produced by PEKA into one table.
    Columns are filenames and rows are k-mers.
    """
    df = pd.DataFrame()
    for f in file_list:
        name = f.split('/')[-1].replace('.tsv', '')
        dft = pd.read_csv(f, sep='\t', index_col=0)
        df[name] = dft['PEKA-score']
    return df


def GetHlines(df_out):
    cumulativeSum = df_out.groupby('cluster').count().cumsum()
    return cumulativeSum.kmers.tolist()


def cluster_kmers(kmers, cluster_on='seq', n_clusters='auto', peka_files=None, results_folder='.', min_token_len=1, weblogos=True, consensus_length='auto'):
    """
    This function clusters k-mers based on their sequence,
    or based on k-mer sequence together with enrichment in PEKA.
    ___________________________________________________________________
    This script takes as inputs:
    1. kmers: a list of k-mers to cluster
    2. an option how to cluster: valid options are 'seq' and 'seq_enrichment'
    3. n_clusters: number of clusters to split k-mers into, can be integer of "auto"
    4. a list of PEKA tsv files, if 'seq_enrichment' is chosen as a mode of
       clustering.

    The script outputs:
    1. a table with grouped k-mers
    2. sequence logo for each k-mer group
    3. distance matrix on which clustering was conducted - this can be used to perform different
      types of clustering outside of this script
    4. a heatmap of clustered k-mers if peka_files are provided
    5. weblogos of motif groups (optional)
    """

    sns.set_style("white")

    # Sort k-mers and remove potential duplicates
    kmers = [k.upper().replace('T', 'U') for k in sorted(set(kmers))]

    # Compute distances between k-mers based on their sequence
    SimilarityMatrixMotifs = clust.GetKmerSimilarityMatrixNoOccurrence(kmers, minLen=min_token_len)
    distanceMatrixMotifs = 1 - SimilarityMatrixMotifs

    # If PEKA files are available, compute zScores
    if peka_files != None:
        PekaScores = GetPEKAScoresFromTsv(peka_files)
        PekaScores.index.name = 'kmers'
        zScores = PekaScores.fillna(0).apply(zscore)

    # Cluster based on sequence
    if cluster_on == 'seq':
        df_out = clust.GetClusteredDf(kmers, similarity_matrix=SimilarityMatrixMotifs, n_clusters=n_clusters)
        dist_matx = pd.DataFrame(distanceMatrixMotifs, index=kmers, columns=kmers)

    # Cluster based on sequence and enrichment in PEKA
    elif cluster_on == 'seq_enrichment':
        # Compute pairwise euclidean distances on PEKA z-scores
        # Perform sorting first
        zScores = zScores.loc[kmers]
        euclideanDistancezScores = sklearn.metrics.pairwise.euclidean_distances(zScores)
        # Compute pairwise k-mer sequence distances
        combined_matrix, weight = clust.CombineDistanceMatrices(distanceMatrixMotifs, euclideanDistancezScores)
        df_out = clust.GetClusteredDf(kmers, similarity_matrix=(1 - combined_matrix), n_clusters=n_clusters)
        dist_matx = pd.DataFrame(combined_matrix, index=kmers, columns=kmers)
    else:
        print('Invalid argument for cluster_on. Must be one of the following: "seq", "seq_enrichment"')

    # Plot heatmap if PEKA files are available
    if peka_files != None:
        df_heatmap = zScores.loc[df_out.kmers.tolist()]
        df_heatmap.to_csv(f'{results_folder}/Heatmap_KmerEnrichmentZScores_NClusters_{str(n_clusters)}_Mode_{str(cluster_on)}.tsv', sep='\t')

        fig = sns.clustermap(df_heatmap, col_cluster=True, row_cluster=False)
        ax = fig.ax_heatmap
        for l in GetHlines(df_out)[:-1]:
            ax.axhline(l, lw=3, color='white')
        fig.savefig(f'{results_folder}/Heatmap_KmerEnrichmentZScores_NClusters_{str(n_clusters)}_Mode_{str(cluster_on)}.pdf', bbox_inches='tight')


    # Add consensus to k-mer clusters
    if consensus_length == 'auto':
        cl = len(kmers[0]) - 1
    else:
        if consensus_length <= len(kmers[0]):
            cl = consensus_length
        else:
            print(f'Warning: Length of consensus is greater than k-mer length. Setting consensus length to k-mer length - {len(kmers[0])} nt')
            cl = len(kmers[0])

    for g, df in df_out.groupby('cluster'):
        motifs = df.kmers.tolist()
        if len(motifs) > 1:
            consensus = align.get_consensus_from_motifs(motifs, l=cl)[0]
        else:
            consensus = motifs[0]
        df_out.loc[df_out.kmers.isin(motifs), 'consensus'] = consensus
        # Plot weblogos if the option is set to True
        if weblogos:
            os.makedirs(f'{results_folder}/weblogos', exist_ok=True)
            if len(motifs) > 1:
                fig, ax = align.PlotKmerLogoLogomaker(motifs, missm_s=-1)
                ax.set_title(consensus)
                sns.despine(right=True, top=True)
                for patch in ax.patches:
                    patch.set_zorder(10)
                ax.tick_params(axis="x", labelbottom=False)
                ax.text(0.8, 0.8, f'n = {len(motifs)}', transform=ax.transAxes)
                fig.savefig(f'{results_folder}/weblogos/Logo_cluster_{g}_{consensus}.pdf', bbox_inches='tight')
    df_out.to_csv(f'{results_folder}/ClusteredKmers_NClusters_{str(n_clusters)}_Mode_{str(cluster_on)}.tsv', sep='\t', index=False)
    dist_matx.to_csv(f'{results_folder}/KmerDistanceMatrix_Mode_{str(cluster_on)}.tsv', sep='\t', header=True, index=True)
    return

def main():
    args = cli()
    cluster_kmers(
    *args
    )

