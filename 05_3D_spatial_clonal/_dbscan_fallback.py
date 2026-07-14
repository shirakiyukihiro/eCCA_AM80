# ============================================================
# DBSCAN implementation (scipy cKDTree; used when scikit-learn is unavailable)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Dependency of the clustering scripts.
# Imported by analyze_tumor_clones_DoG_Zcrop.py and test_deconvolution_deep_zone.py
# when scikit-learn is unavailable.
# ============================================================

"""Self-contained DBSCAN (scipy.cKDTree only; scikit-learn not required)."""
import numpy as np
from scipy.spatial import cKDTree

def dbscan_kdtree(points, eps, min_samples):
    """
    DBSCAN clustering.

    points      : (N, 3) array of coordinates
    eps         : neighbourhood radius
    min_samples : minimum number of neighbours (including the point itself)
                  for a point to be a core point

    Returns labels (N,); -1 = noise, 0+ = cluster id.
    """
    n = len(points)
    if n == 0:
        return np.array([], dtype=int)
    tree = cKDTree(points)
    neighbors = tree.query_ball_point(points, r=eps)
    is_core = np.array([len(ns) >= min_samples for ns in neighbors])
    labels = np.full(n, -1, dtype=int)
    cluster_id = 0
    visited = np.zeros(n, dtype=bool)
    for i in range(n):
        if visited[i] or not is_core[i]:
            continue
        # BFS
        queue = [i]
        while queue:
            j = queue.pop()
            if visited[j]:
                continue
            visited[j] = True
            labels[j] = cluster_id
            if is_core[j]:
                for k in neighbors[j]:
                    if not visited[k]:
                        queue.append(k)
        cluster_id += 1
    return labels
