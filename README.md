# cluster_kmers
Get clusters of k-mers based solely on their sequence or in combination with enrichment in PEKA. 

## Installation

```pip install git+https://github.com/ulelab/cluster_kmers.git@master```


## Usage
```
cluster_kmers [-h] -k KMERS [KMERS ...] -o RESULTS_FOLDER [-co {seq,seq_enrichment}] [-n N_CLUSTERS] [-peka PEKA_FILES [PEKA_FILES ...]] [-tl MIN_TOKEN_LENGTH] [-wl {True,False}] [-cons CONSENSUS_LENGTH]
```

Cluster k-mers based on sequence or on a combination of sequence and enrichment in CLIP data.
```
required arguments:
  -k KMERS [KMERS ...], --kmers KMERS [KMERS ...]
                        A list of k-mers encoded in RNA alphabet separated by spaces: AAGG GGAG GCCU.
  -o RESULTS_FOLDER, --output_folder RESULTS_FOLDER
                        A path to an existing output folder, for example "~/results"

optional arguments:
  -co {seq,seq_enrichment}, --cluster_on {seq,seq_enrichment}
                        Inputs to clustering. Valid options are: seq - Cluster only based on sequence similarity. seq_enrichment - Cluster based on sequence similarity and based on enrichment of motifs in CLIP data (this option requires arguments passed to -peka)
  -n N_CLUSTERS, --n_clusters N_CLUSTERS
                        Number of clusters to split k-mers into. Valid options are "auto" or integer. Default is "auto".
  -peka PEKA_FILES [PEKA_FILES ...], --peka_files PEKA_FILES [PEKA_FILES ...]
                        A list of peka output files with extensions *mer_distribution_{region_name}.tsv, separated by spaces.
  -tl MIN_TOKEN_LENGTH, --min_token_length MIN_TOKEN_LENGTH
                        Minimal length of a substrings used for clustering. For k-mers with lengths greater than 5, setting this value to be greater than 1 can improve the results of clustering.
  -wl {True,False}, --weblogos {True,False}
                        Whether to plot weblogos for motif groups, True by default.
  -cons CONSENSUS_LENGTH, --consensus_length CONSENSUS_LENGTH
                        Length of consensus sequence to name k-mer groups. Automatically this length is determined as k-mer length - 1. Valid choices are "auto" or integer.
```
