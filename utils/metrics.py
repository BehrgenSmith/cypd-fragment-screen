import numpy as np
from typing import Sequence, Dict
from scipy.spatial.distance import pdist
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage,fcluster

def normalized_auc(x, y):
    """
    Compute the AUC of a curve using the trapezoidal rule,
    normalized so the result is between 0 and 1 (unit square).

    Parameters
    ----------
    x : array-like
        X-values of the curve (must be monotonic).
    y : array-like
        Y-values of the curve.

    Returns
    -------
    float
        Normalized AUC (0 to 1).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Sort by x
    idx = np.argsort(x)
    x, y = x[idx], y[idx]

    # Normalize to [0, 1]
    x_norm = (x - x.min()) / (x.max() - x.min())
    y_norm = (y - y.min()) / (y.max() - y.min())

    # Compute AUC on normalized data
    auc = np.trapezoid(y_norm, x_norm)

    return auc



def rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Compute the root-mean-square deviation (RMSD) between two coordinate sets.

    Parameters
    ----------
    coords1, coords2 : np.ndarray
        Arrays of shape (N, 3) representing 3D coordinates.
        They must correspond atom-by-atom (same order, same size).

    Returns
    -------
    float
        RMSD value in Å.
    """
    if coords1.shape != coords2.shape:
        raise ValueError(f"Coordinate shapes differ: {coords1.shape} vs {coords2.shape}")
    diff = coords1 - coords2
    return np.sqrt((diff ** 2).sum() / coords1.shape[0])

def pairwise_rmsd(arrays_a, arrays_b):
    """
    Compute pairwise RMSD between two lists/arrays of conformations.

    Parameters
    ----------
    arrays_a : list[np.ndarray] or np.ndarray
        Shape (M, N, 3) or list of M arrays (N, 3).
    arrays_b : list[np.ndarray] or np.ndarray
        Shape (K, N, 3) or list of K arrays (N, 3).

    Returns
    -------
    np.ndarray
        RMSD matrix of shape (M, K)
    """
    arrays_a = np.array(arrays_a)
    arrays_b = np.array(arrays_b)
    m, k = len(arrays_a), len(arrays_b)

    rmsd_mat = np.zeros((m, k))
    for i in range(m):
        for j in range(k):
            rmsd_mat[i, j] = rmsd(arrays_a[i], arrays_b[j])
    return rmsd_mat

def cluster_counts_over_similarity(
    fps: Sequence[Sequence[int]],
    n_divisions: int = 10,
    metric: str = "jaccard",
    method: str = "average",
    return_linkage = False
) -> Dict[int, int]:
    """
    For a set of fingerprints, compute the number of clusters at a series of
    similarity thresholds in (0, 1].

    Returns a dict mapping integer similarity * 100 -> n_clusters.
    e.g. similarity 0.1 -> key 10, similarity 0.25 -> key 25.
    """
    fps_arr = np.asarray(fps, dtype=bool)

    # pairwise distances + linkage (done once)
    dists = pdist(fps_arr, metric=metric)
    Z = linkage(dists, method=method, optimal_ordering=True)

    # similarity thresholds in (0,1], e.g. 0.1,0.2,... for n_divisions=10
    sims = np.linspace(0, 1, n_divisions + 1)[1:]
    cluster_counts = []
    for s in tqdm(sims):
        d = 1.0 - s  # distance cutoff corresponding to similarity s
        labels = fcluster(Z, t=d, criterion="distance")
        n_clusters = np.unique(labels).size
        key = int(round(s * 100))  # e.g. 0.1 -> 10, 0.25 -> 25
        cluster_counts.append(n_clusters)


    if return_linkage:
        np.array(cluster_counts), sims, Z, dists
    else:

        return np.array(cluster_counts), sims