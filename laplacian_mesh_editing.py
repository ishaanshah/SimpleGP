import numpy as np
import networkx as nx
import time
import polyscope as ps
from trimesh import exchange
from trimesh import graph
from argparse import ArgumentParser
from scipy import sparse

parser = ArgumentParser()
parser.add_argument("--file", default="./spot/spot_quadrangulated.obj", type=str)
parser.add_argument("--boundary_file", default="./head_boundary", type=str)
parser.add_argument("--handle_idx", default=63, type=int)
parser.add_argument("--handle_delta", default=[0, 0.1, 0], nargs=3)

args = parser.parse_args()

ps.init()
mesh = exchange.load.load(args.file)

# Create graph from mesh
G = graph.vertex_adjacency_graph(mesh)

# Define boundary
boundary_verts = []
with open(args.boundary_file, "r") as f:
    for line in f.readlines():
        boundary_verts.append(int(line))

boundary_verts = np.array(boundary_verts)
boundary = []
for i in range(len(boundary_verts)):
    src = boundary_verts[i]
    trg = boundary_verts[(i+1) % len(boundary_verts)]
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
src_mesh.add_color_quantity("edit_params", colors, defined_on="vertices", enabled=True)

# Calculate laplacian coordinates
L = nx.laplacian_matrix(G_sub)
D_inv = sparse.diags(1 / np.asarray(nx.adjacency_matrix(G_sub).sum(axis=0)).flatten())

# Find laplacian coordinates
G_sub_nodes = np.asarray(G_sub.nodes)
pos = mesh.vertices[G_sub_nodes]
delta = L @ pos

# Plot delta
plot_delta = np.zeros((len(mesh.vertices), 3))
plot_delta[G_sub_nodes] = delta
src_mesh.add_vector_quantity("delta", plot_delta, defined_on="vertices", enabled=True, vectortype="ambient", radius=0.001, color=[1, 0, 1])

ps.show()
