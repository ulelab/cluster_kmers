import pandas as pd
import numpy as np
import textdistance as td
import sklearn.cluster


def get_unique_string_tokensNoOccurrence(test_str):
    tokens = [test_str[i: j] for i in range(len(test_str)) for j in range(i + 1, len(test_str) + 1)]
    # Filter tokens that correspond to test_str
    tokens = [t for t in tokens if t != test_str]
    unique_tokens = sorted(list(set(tokens)))
    return unique_tokens


def SortTokensByLength(tokens):
    lenDict = {}
    for length in range(np.min([len(t) for t in tokens]), np.max([len(t) for t in tokens]) + 1):
        lenDict[length-1] = [t for t in tokens if len(t) == length]
    return lenDict


def ComputePairwiseKmerSimilarityNoOccurrence(k1, k2, minLen=1):
    dict1 = SortTokensByLength(get_unique_string_tokensNoOccurrence(k1))
    dict2 = SortTokensByLength(get_unique_string_tokensNoOccurrence(k2))
    jaccardSimilarities = []
    lengths = [l for l in dict1.keys() if l in dict2.keys()]
    for tokenLength in lengths:
        if tokenLength >= minLen:
            jaccardSimilarities.append(td.jaccard.similarity(dict1[tokenLength], dict2[tokenLength]))
    # Mean Jaaccard is between 0 and 1
    MeanJaccard = np.mean(jaccardSimilarities)
    return MeanJaccard


def GetKmerSimilarityMatrixNoOccurrence(words, minLen=1):
    words = np.asarray(words)
    SimilarityMatrix = np.array(
        [[ComputePairwiseKmerSimilarityNoOccurrence(w1, w2, minLen) for w1 in words] for w2 in words]
        )
    return SimilarityMatrix

def CombineDistanceMatrices(m1, m2, weight=1):
    """
    Combines two distances matrices for orthogonal variables.
    """
    # Standard scaling over columns
    scaler = sklearn.preprocessing.StandardScaler()
    m1Scaled = scaler.fit_transform(m1)
    m2Scaled = scaler.fit_transform(m2)
    # Apply min max scaling to move the matrices to be between 0 and 1
    MinMaxScaler = sklearn.preprocessing.MinMaxScaler()
    m1Normalized = MinMaxScaler.fit_transform(m1Scaled)
    m2Normalized = MinMaxScaler.fit_transform(m2Scaled)
    # Combine into one matrix using Pythagoras
    # Proportional contribution of m1 on the scale of 0 to 1
    m1Weighted = m1Normalized * weight
    m2Weighted = m2Normalized * 1
    combined_matrix = np.hypot(m1Weighted, m2Weighted)
    return combined_matrix, weight


def GetClusteredDf(kmers, similarity_matrix, n_clusters):
    # Cluster with affinity propagation if n_clusters='auto'
    if n_clusters == 'auto':
        affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=0.5, max_iter=1000, convergence_iter=200, random_state=0)
        X = affprop.fit(similarity_matrix)
        df = pd.DataFrame([kmers, X.labels_], index=['kmers', 'cluster']).T.sort_values(by='cluster')
     # Cluster with KMeans if n_clusters is an integer
    elif type(n_clusters) == int:
        if n_clusters > len(kmers):
            print(f'Number of specified clusters is greater than the number of k-mers passed ({len(kmers)}).')
            print('Exiting.')
            pass
        else:
            distanceMatrix = 1 - similarity_matrix
            kmeans = sklearn.cluster.KMeans(n_clusters=n_clusters,  init='k-means++', n_init=10)
            X = kmeans.fit(distanceMatrix)
            df = pd.DataFrame([kmers, X.labels_], index=['kmers', 'cluster']).T.sort_values(by='cluster')
    else:
        print('Invalid argument for n_clusters. Must be one of the following: "auto" or integer')
    return df


