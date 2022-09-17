import numpy as np
import polyscope as ps
import time
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from trimesh import exchange
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--file", default="./spot/pcd_1k.xyz")
parser.add_argument("--kd_tree", action="store_true")
parser.add_argument("--algorithm", default="svd", choices=["svd", "lsqr_point", "lsqr_plane"])
parser.add_argument("--iters", default=500, type=int)

args = parser.parse_args()

ps.init()
spot = exchange.load.load(args.file)

src = np.asarray(spot.vertices)
src_pcd = ps.register_point_cloud("source", src)

rot = Rotation.from_euler("XYZ", [30, 0, 0], degrees=True).as_matrix()

trg = src @ rot + np.random.rand()
trg_pcd = ps.register_point_cloud("target", trg)

t_m = trg.mean(axis=0)
t_mc = trg - t_m
kd_tree = KDTree(t_mc)
t0 = time.time()
for i in range(args.iters):
    s_m = src.mean(axis=0)
    s_mc = src - s_m
    if args.kd_tree:
        _, corr_idx = kd_tree.query(s_mc)
    else:
        corr_idx = np.argmin(cdist(s_mc, t_mc), axis=1)

    cov = t_mc[corr_idx].T @ s_mc
    u, _, vt = np.linalg.svd(cov)
    r = u @ vt
    if not np.isclose(np.linalg.det(r), 1):
        r = u @ np.diag([1, 1, -1]) @ vt

    t = t_m - (rot @ s_m.T).T
    src = (rot @ src.T).T + t
    if np.isclose(src, trg).all():
        break
    print(np.mean(np.linalg.norm(src - trg)))
    
pred_pcd = ps.register_point_cloud("predicted", src)
print("Time take:", time.time() - t0)

ps.show()
