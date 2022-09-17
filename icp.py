import numpy as np
import polyscope as ps
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from trimesh import exchange

ps.init()
spot = exchange.load.load("./spot/pcd_1k.xyz")

src = np.asarray(spot.vertices)
src_pcd = ps.register_point_cloud("source", src)

rot = Rotation.from_euler("XYZ", [5, 0, 0], degrees=True).as_matrix()

trg = src @ rot + np.ones(3)
trg_pcd = ps.register_point_cloud("target", trg)

for i in range(100):
    s_m, t_m = src.mean(axis=0), trg.mean(axis=0)
    s_mc, t_mc = src - s_m, trg - t_m
    corr_idx = np.argmin(cdist(s_mc, t_mc), axis=1)

    cov = t_mc[corr_idx].T @ s_mc
    u, _, vt = np.linalg.svd(cov)
    r = u @ vt
    if not np.isclose(np.linalg.det(r), 1):
        r = u @ np.diag([1, 1, -1]) @ vt

    t = t_m - s_m @ rot
    src = src @ rot + t
    if np.linalg.norm(src - trg) < 1e-5:
        break
    print(np.mean(np.linalg.norm(src - trg)))
    
print(np.isclose(src, trg).all())
pred_pcd = ps.register_point_cloud("predicted", src)

ps.show()
