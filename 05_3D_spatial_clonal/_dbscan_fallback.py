# ============================================================
# DBSCAN implementation (scipy cKDTree; used when scikit-learn is unavailable)
# Study      : Distinct CAF lineages and an Am80 stromal-reprogramming strategy
#              in extrahepatic cholangiocarcinoma (Meflin/ISLR).
# Repository : eCCA_AM80  (https://github.com/shirakiyukihiro/eCCA_AM80)
# Dependency of the clustering scripts.
# Paths made relative for public release: raw .ims files in ./data (not distributed;
# available on request), region masks in ./data/mask_polygons, results in ./output.
# ============================================================

"""DBSCAN 自前実装 (scipy.cKDTree 使用、sklearn 不要)"""
import numpy as np
from scipy.spatial import cKDTree

def dbscan_kdtree(points, eps, min_samples):
    """
    DBSCAN 同等の clustering
    points: (N, 3) array
    eps: 近傍半径
    min_samples: core point の最小近傍数 (自分含む)

    返り値: labels (N,), -1 = noise, 0+ = cluster id
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
