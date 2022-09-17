import numpy as np
import polyscope as ps
import time
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from scipy.optimize import least_squares
from trimesh import exchange
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--file", default="./spot/pcd_1k.xyz")
parser.add_argument("--kd_tree", action="store_true")
parser.add_argument("--algorithm", default="svd", choices=["svd", "lsqr_point", "lsqr_plane"])
parser.add_argument("--iters", default=100, type=int)

args = parser.parse_args()

def svd_icp(s_mc, t_mc, s_m, t_m, corr_idx):
    cov = t_mc[corr_idx].T @ s_mc
    u, _, vt = np.linalg.svd(cov)
    r = u @ vt
    if not np.isclose(np.linalg.det(r), 1):
        r = u @ np.diag([1, 1, -1]) @ vt

    t = t_m - (r @ s_m.T).T
    return r, t

def lsqr_icp(src, trg, corr_idx):
    def point_fwd(x):
        r = Rotation.from_quat(x[:4]).as_matrix()
        t = x[4:]
        return np.sum(((r @ src.T).T + t - trg[corr_idx]) ** 2, axis=1)

    if args.algorithm == "lsqr_point":
        x = least_squares(point_fwd, np.random.uniform(size=7))
    else:
        x = least_squares(point_fwd, np.random.uniform(size=7))

    r = Rotation.from_quat(x.x[:4]).as_matrix()
    t = x.x[4:]

    return r, t

ps.init()
spot = exchange.load.load(args.file)

src = np.asarray(spot.vertices)
src_pcd = ps.register_point_cloud("source", src)

rot = Rotation.from_euler("XYZ", [30, 60, 0], degrees=True).as_matrix()

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

    if args.algorithm == "svd":
        r, t = svd_icp(s_mc, t_mc, s_m, t_m, corr_idx)
    else:
        r, t = lsqr_icp(src, trg, corr_idx)

    src = (r @ src.T).T + t
    if np.isclose(src, trg).all():
        break
    print(np.mean(np.linalg.norm(src - trg)))
    
pred_pcd = ps.register_point_cloud("predicted", src)
print("Time taken:", time.time() - t0)

ps.show()
