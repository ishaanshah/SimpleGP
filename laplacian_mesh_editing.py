import numpy as np
import networkx as nx
import time
import polyscope as ps
from trimesh import exchange
from trimesh import graph
from argparse import ArgumentParser
from scipy import sparse

# Debug variables
debug_ring = False

parser = ArgumentParser()
parser.add_argument("--file", default="./spot/spot_quadrangulated.obj", type=str)
parser.add_argument("--boundary_file", default="./head_boundary", type=str)
parser.add_argument("--handle_idx", default=1453, type=int)
parser.add_argument("--handle_delta", default=[0.2, 0.1, 0.1], nargs=3)

args = parser.parse_args()

ps.init()
mesh = exchange.load.load(args.file)

# Create graph from mesh
G = graph.vertex_adjacency_graph(mesh)
G_pos = mesh.vertices

# Define boundary
boundary_verts = []
with open(args.boundary_file, "r") as f:
    for line in f.readlines():
        boundary_verts.append(int(line))

boundary_verts = np.array(boundary_verts)
boundary = []
for i in range(len(boundary_verts)):
    src = boundary_verts[i]
    trg = boundary_verts[(i + 1) % len(boundary_verts)]
    path = nx.shortest_path(G, source=src, target=trg)
    boundary += path[:-1]
boundary = np.asarray(boundary)

# Register source mesh
src_mesh = ps.register_surface_mesh("source", mesh.vertices, mesh.faces)
src_mesh.set_smooth_shade(True)

# Get editable subgraph
G_sub = G.copy()
G_sub.remove_nodes_from(boundary)
editable = np.array(list(nx.node_connected_component(G_sub, args.handle_idx)))
G_sub = G.subgraph(np.concatenate((boundary, editable)))

# Mark boundary, handle and editable vertices
colors = np.zeros((len(mesh.vertices), 3))
colors[editable] = [0, 1, 0]
colors[boundary] = [1, 0, 0]
colors[args.handle_idx] = [0, 0, 1]
src_mesh.add_color_quantity("edit_params", colors, defined_on="vertices", enabled=False)

# Calculate laplacian coordinates
L = nx.laplacian_matrix(G_sub)
D_inv = sparse.diags(1 / np.asarray(nx.adjacency_matrix(G_sub).sum(axis=0)).flatten())
L = D_inv @ L

# Find laplacian coordinates
G_sub_nodes = np.asarray(G_sub.nodes)
G_sub_cnt = G_sub_nodes.shape[0]
G_sub_pos = mesh.vertices[G_sub_nodes]
delta = L @ G_sub_pos
L = L.todense()

# Get mapping from node_id to idx
id2idx = {}
for idx, id in enumerate(G_sub_nodes):
    id2idx[id] = idx

# Plot delta
plot_delta = np.zeros((len(mesh.vertices), 3))
plot_delta[G_sub_nodes] = delta
src_mesh.add_vector_quantity(
    "delta",
    plot_delta,
    defined_on="vertices",
    enabled=False,
    vectortype="ambient",
    radius=0.001,
    color=[1, 0, 1],
)

# Build linear system to solve
system = np.zeros((3 * G_sub_cnt, 3 * G_sub_cnt))
system[0:G_sub_cnt, 0:G_sub_cnt] = -L
system[G_sub_cnt : 2 * G_sub_cnt, G_sub_cnt : 2 * G_sub_cnt] = -L
system[2 * G_sub_cnt : 3 * G_sub_cnt, 2 * G_sub_cnt : 3 * G_sub_cnt] = -L

for i in range(G_sub_cnt):
    node = G_sub_nodes[i]

    # Get neighbor in mesh idx
    neigh = np.asarray(list(G_sub.neighbors(node)))

    # Get idx of neighbour in subgraph
    n_idx = np.zeros_like(neigh)
    for j in range(neigh.shape[0]):
        n_idx[j] = id2idx[neigh[j]]

    # Debug
    if debug_ring:
        ring = np.zeros((len(mesh.vertices), 3))
        ring[node] = [0, 0, 1]
        ring[neigh] = [1, 0, 0]
        src_mesh.add_color_quantity("ring", ring, defined_on="vertices", enabled=True)
        break

    n_idx = np.concatenate((np.asarray([i]), n_idx))
    A = np.zeros([n_idx.shape[0] * 3, 7])
    for idx, j in enumerate(n_idx):
        A[idx] = [G_sub_pos[j, 0], 0, G_sub_pos[j, 2], -G_sub_pos[j, 1], 1, 0, 0]
        A[idx + n_idx.shape[0]] = [
            G_sub_pos[j, 1],
            -G_sub_pos[j, 2],
            0,
            G_sub_pos[j, 0],
            0,
            1,
            0,
        ]
        A[idx + 2 * n_idx.shape[0]] = [
            G_sub_pos[j, 2],
            G_sub_pos[j, 1],
            -G_sub_pos[j, 0],
            0,
            0,
            0,
            1,
        ]

    Ai = np.linalg.pinv(A)
    s = Ai[0]
    h = Ai[1:4]
    t = Ai[4:]

    T = np.zeros((3, 3, s.shape[0]))
    T[0, 0] = s
    T[1, 1] = s
    T[2, 2] = s
    T[0, 1] = -h[2]
    T[0, 2] = h[1]
    T[1, 0] = h[2]
    T[1, 2] = -h[0]
    T[2, 0] = -h[1]
    T[2, 1] = h[0]
    T_delta = T.transpose(0, 2, 1) @ delta[i].T

    system[
        i, np.concatenate((n_idx, n_idx + G_sub_cnt, n_idx + 2 * G_sub_cnt))
    ] += T_delta[0]
    system[
        i + G_sub_cnt, np.concatenate((n_idx, n_idx + G_sub_cnt, n_idx + 2 * G_sub_cnt))
    ] += T_delta[1]
    system[
        i + 2 * G_sub_cnt,
        np.concatenate((n_idx, n_idx + G_sub_cnt, n_idx + 2 * G_sub_cnt)),
    ] += T_delta[2]

b = np.zeros(3 * G_sub_cnt)

# # Add constraints boundary
b_idx = np.zeros_like(boundary_verts)
for i in range(boundary_verts.shape[0]):
    b_idx[i] = id2idx[boundary_verts[i]]
for i in b_idx:
    Ln = np.zeros((3, 3 * G_sub_cnt))
    Ln[0, i] = 1
    Ln[1, i + G_sub_cnt] = 1
    Ln[2, i + 2 * G_sub_cnt] = 1
    system = np.vstack((system, Ln))
    b = np.concatenate((b, G_sub_pos[i]))

# Add constraints for handle
h_idx = id2idx[args.handle_idx]
Ln = np.zeros((3, 3 * G_sub_cnt))
Ln[0, h_idx] = 1
Ln[1, h_idx + G_sub_cnt] = 1
Ln[2, h_idx + 2 * G_sub_cnt] = 1

A = np.vstack((system, Ln))
b = np.concatenate((b, G_sub_pos[h_idx] + np.array(args.handle_delta)))

# Solve equation Ax = b (x = v')
x = sparse.linalg.lsqr(sparse.coo_matrix(A), b)[0]

# Update mesh positions
vp = x.reshape(3, -1).T
n_mesh = mesh.copy()
n_mesh.vertices[G_sub_nodes] = vp
n_mesh.vertices += np.asarray([1, 0, 0])
new_mesh = ps.register_surface_mesh(
    "edited", n_mesh.vertices, n_mesh.faces, color=[0.5, 0, 0]
)
new_mesh.set_smooth_shade(True)

ps.show()
