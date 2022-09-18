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
parser.add_argument(
    "--algorithm", default="svd", choices=["svd", "lsqr_point", "lsqr_plane"]
)
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


def lsqr_icp(src, trg, corr_idx, trg_normals):
    def point_fwd(x):
        r = Rotation.from_quat(x[:4]).as_matrix()
        t = x[4:]
        return np.sum(((r @ src.T).T + t - trg[corr_idx]) ** 2, axis=1)

    def plane_fwd(x):
        r = Rotation.from_quat(x[:4]).as_matrix()
        t = x[4:]
        trg_n = np.array(trg_normals)[corr_idx]
        return np.sum((((r @ src.T).T + t) - trg[corr_idx]) * trg_n, axis=1)

    if args.algorithm == "lsqr_point":
        x = least_squares(point_fwd, np.random.uniform(size=7))
    else:
        x = least_squares(plane_fwd, np.random.uniform(size=7))

    r = Rotation.from_quat(x.x[:4]).as_matrix()
    t = x.x[4:]

    return r, t


ps.init()
pcd = exchange.load.load(args.file)

# Load normals (trimesh doesn't support loading normals from ".xyz")
# This makes the program only work with ".xyz" files, however this can
# be easily fixed by estimating the normals
src_normals = []
with open(args.file) as f:
    for line in f.readlines():
        src_normals.append(line.split()[3:])

src_normals = np.array(src_normals, dtype=float)

src = np.asarray(pcd.vertices)
src_pcd = ps.register_point_cloud("source", src)
src_pcd.add_vector_quantity("normals", src_normals)

rot = Rotation.from_euler("XYZ", np.random.rand(3) * 90, degrees=True).as_matrix()

trg = src @ rot + np.random.rand()
trg_normals = src_normals @ rot
trg_pcd = ps.register_point_cloud("target", trg)
trg_pcd.add_vector_quantity("normals", trg_normals)

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
        r, t = lsqr_icp(src, trg, corr_idx, trg_normals)

    src = (r @ src.T).T + t
    if np.isclose(src, trg).all():
        break
    print(np.mean(np.linalg.norm(src - trg)))

pred_pcd = ps.register_point_cloud("predicted", src)
print("Time taken:", time.time() - t0)

ps.show()
