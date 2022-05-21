import copy
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
from casm.info import get_supercell_info
from casm_figures.cluster import select_cluster_sites
from casm_figures.lattice import plot_lattice_cell, \
    plot_cartesian_axes, apply_cabinet
from casm_figures.view import make_cartesian_view_basis, make_view_basis
from casm_figures.supercell import supercell_sites_df

plt.rcParams.update({'font.family': 'Times New Roman'})
plt.rcParams.update({'font.size': 10})
plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})

root = "hex_binary_2d"
clust_json_path = "hex_binary_2d/basis_sets/bset.default/clust.json"
x_stop = 30 # larger -> wider figure
figure_inches_wide = 16
figure_inches_tall = 12

dim = 2

# Read prim, prim lattice
prim = None
L = None    # prim lattice column matrix
with open(root + '/.casm/prim.json','r') as f:
    prim = json.load(f)
    L = np.array(prim["lattice_vectors"]).transpose()

# Get prim cart coordinates
I = np.eye(3, dtype=int)
info = get_supercell_info(prim=prim,
                          transformation_matrix_to_super=I.tolist(),
                          properties=["cart_coordinate"])
prim_cart_basis = np.array(info["cart_coordinate"]).transpose()


# Set plotting parameters for sites by asymmetric unit
s_scale = 40
asym_scatter_kwargs = {
    0: {
        'c': 'red',
        's': 1.*s_scale,
        'edgecolors': 'black'
    },
    1: {
        'c': 'blue',
        's': 1.*s_scale,
        'edgecolors': 'black'
    },
    2: {
        'c': 'white',
        's': 2.*s_scale,
        'edgecolors': 'black'
    }
}

def view_sites(view_basis, _df, asym_scatter_kwargs={}, cabinet=None):
    df = _df.copy()
    view_basis_inv = np.linalg.pinv(view_basis)
    values = view_basis_inv @ df[['x','y','z']].values.transpose()
    apply_cabinet(cabinet, values)
    df['px'] = values[0,:]
    df['py'] = values[1,:]
    df['pz'] = values[2,:]
    df = df.sort_values('pz')
    df['c'] = [asym_scatter_kwargs[asym]['c'] for asym in df['asym']]
    df['s'] = [asym_scatter_kwargs[asym]['s'] for asym in df['asym']]
    df['edgecolors'] = [asym_scatter_kwargs[asym]['edgecolors'] for asym in df['asym']]
    return df

def plot_background(view_basis, L, df, T, supercell_shift_ijk, cabinet=None, T_conventional=None, axis_labels=['a', 'b', 'c'], hex=False):
    # Plot cell
    unitcell_cell = True
    supercell_cell = False
    cell_plot_kwargs = {
        'color': 'black',
        'linestyle': 'dashed'}
    cell_arrows = True
    cell_arrow_kwargs = {
        'head_width':.1,
        'length_includes_head':True,
        'color':'black'
    }
    cell_annotate = True
    cell_annotate_kwargs = {
        'textcoords': 'offset points',
        'xytext': (-10, -10)
        }
    if T_conventional is None:
        T_conventional = np.eye(3, dtype=int)

    L_unitcell = None
    unitcell_shift = None
    supercell_shift = None
    if hex:
        L_unitcell = L
        hex_shift = np.array([1, 1, 0], dtype=int)
        unitcell_shift = L @ (supercell_shift_ijk + hex_shift)
        supercell_shift = L @ supercell_shift_ijk
    else:
        L_unitcell = L @ T_conventional
        unitcell_shift = L @ supercell_shift_ijk
        supercell_shift = L @ supercell_shift_ijk

    # plot unit cell
    plot_lattice_cell(view_basis, L_unitcell, cabinet=cabinet,
                      cell=unitcell_cell,
                      plot_kwargs=cell_plot_kwargs,
                      arrows=cell_arrows,
                      arrow_kwargs=cell_arrow_kwargs,
                      annotate=cell_annotate,
                      annotate_kwargs=cell_annotate_kwargs, dim=dim,
                      axis_labels=axis_labels,
                      shift=unitcell_shift,
                      hex=hex)
    # plot super cell
    plot_lattice_cell(view_basis, L @ T, cabinet=cabinet,
                      cell=supercell_cell,
                      plot_kwargs=cell_plot_kwargs,
                      arrows=False,
                      annotate=False, dim=dim,
                      shift=supercell_shift)

    # Plot sites
    N = len(df['px'])
    for i in range(N):
        plt.scatter(x=[df['px'].values[i]], y=[df['py'].values[i]], c=[df['c'].values[i]],
            s=[df['s'].values[i]], edgecolors=df['edgecolors'].values[i], linewidths=1,
            zorder=df['pz'].values[i])

    # Plot Cartesian axes
    cartesian_axes_scale = 1.0
    cartesian_axes_head_width = cartesian_axes_scale / 10.
    cartesian_axes_scale_origin = L @ (supercell_shift_ijk - np.array([3, 0, 1]))
    cartesian_axes_arrow_kwargs = {
        'color': 'black',
        'head_width':cartesian_axes_head_width,
        'length_includes_head':True}
    plot_cartesian_axes(view_basis, scale=cartesian_axes_scale, cabinet=cabinet,
        shift=cartesian_axes_scale_origin,
        arrow_kwargs=cartesian_axes_arrow_kwargs,
        dim=dim)

    # Axes settings
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_axis_off()

def get_cluster_box_2d(sites):
    if len(sites) == 0:
        return (np.eye(3, dtype=int), np.array([0, 0, 0], dtype=int))
    big = 100
    i_min = big
    i_max = -big
    j_min = big
    j_max = -big
    k_min = big
    k_max = -big
    for site in sites:
        ijk = np.array(site[-3:], dtype=int)

        if ijk[0] < i_min:
            i_min = ijk[0]
        if ijk[0] > i_max:
            i_max = ijk[0]
        if ijk[1] < j_min:
            j_min = ijk[1]
        if ijk[1] > j_max:
            j_max = ijk[1]
        if ijk[2] < k_min:
            k_min = ijk[2]
        if ijk[2] > k_max:
            k_max = ijk[2]

    T = np.eye(3, dtype=int)
    T[0,0] = i_max - i_min + 1
    T[1,1] = j_max - j_min + 1
    T[2,2] = k_max - k_min + 1

    return (T, np.array([i_min, j_min, k_min], dtype=int))

def get_cluster_box(L, prim_cart_basis, sites, T_conventional=None):
    if T_conventional is None:
        T_conventional = np.eye(3, dtype=int)
    if len(sites) == 0:
        return (np.eye(3, dtype=int), np.array([0, 0, 0], dtype=int))
    L_conventional_inv = np.linalg.pinv(L @ T_conventional)
    eps = 1e-5
    big = 100
    i_min = big
    i_max = -big
    j_min = big
    j_max = -big
    k_min = big
    k_max = -big
    for site in sites:
        b = site[0]
        ijk = np.array(site[-3:], dtype=int)
        cart = prim_cart_basis[:,b] + L @ ijk

        # Frac (in conventional cell)
        frac = L_conventional_inv @ cart

        if frac[0] < i_min:
            i_min = frac[0]
        if frac[0] > i_max:
            i_max = frac[0]
        if frac[1] < j_min:
            j_min = frac[1]
        if frac[1] > j_max:
            j_max = frac[1]
        if frac[2] < k_min:
            k_min = frac[2]
        if frac[2] > k_max:
            k_max = frac[2]

    T2 = np.eye(3, dtype=int)
    print(i_min, i_max)
    print(j_min, j_max)
    print(k_min, k_max)
    T2[0,0] = math.ceil(i_max - eps) - math.floor(i_min + eps)
    T2[1,1] = math.ceil(j_max - eps) - math.floor(j_min + eps)
    T2[2,2] = math.ceil(k_max - eps) - math.floor(k_min + eps)
    for i in range(3):
        if T2[i][i] < 1:
            T2[i][i] = 1
    shift2_ijk = np.array([math.floor(i_min + eps), math.floor(j_min + eps), math.floor(k_min + eps)], dtype=int)

    T = T_conventional @ T2
    shift_ijk = T_conventional @ shift2_ijk

    return (T, shift_ijk)

def plot_cluster(view_basis, L, prim_cart_basis, orbit_index, sites,
    cluster_shift_ijk=None, zorder=-1000, cabinet=None,
    cluster_scatter_kwargs=None, cluster_line_kwargs=None,
    cluster_annotate_kwargs=None):

    view_basis_inv = np.linalg.pinv(view_basis)
    eps = 1e-5
    if cluster_annotate_kwargs is None:
        cluster_annotate_kwargs = {
            'textcoords': 'offset points',
            'xytext': (-15, 5),
            'zorder': 1000
            }
    if cluster_line_kwargs is None:
        cluster_line_kwargs = {
            'linewidth': 6
        }
    if cluster_scatter_kwargs is None:
        cluster_scatter_kwargs = {
            's': 4 * s_scale
        }
    if cluster_shift_ijk is None:
        cluster_shift_ijk = np.array([0, 0, 0], dtype=int)
    label_xy = None

    i = 0
    while i < len(sites):
        site_i = sites[i]
        ijk_i = np.array(site_i[-3:], dtype=int) + cluster_shift_ijk
        cart_i = prim_cart_basis[:, site_i[0]] + L @ ijk_i
        values_i = np.zeros((3,1))
        values_i[:,0] = view_basis_inv @ cart_i
        apply_cabinet(cabinet, values_i)

        plt.scatter(x=[values_i[0]], y=[values_i[1]], c=[cluster_color],
            s=[cluster_scatter_kwargs['s']], linewidth=None, zorder=zorder)

        # try to find upper left site in cluster to label
        site_xy = (values_i[0], values_i[1])
        if i == 0:
            label_xy = site_xy
        else:
            if values_i[0] < (label_xy[0] - eps):
                label_xy = site_xy
            elif math.isclose(label_xy[0], values_i[0]) and values_i[1] > label_xy[1] + eps:
                label_xy = site_xy

        j = i + 1
        while j < len(sites):
            site_j = sites[j]
            ijk_j = np.array(site_j[-3:], dtype=int) + cluster_shift_ijk
            cart_j = prim_cart_basis[:, site_j[0]] + L @ ijk_j
            values_j = np.zeros((3,1))
            values_j[:,0] = view_basis_inv @ cart_j
            apply_cabinet(cabinet, values_j)

            kwargs = copy.copy(cluster_line_kwargs)
            kwargs['color'] = cluster_color
            kwargs['zorder'] = zorder

            x = values_i[0]
            y = values_i[1]
            dx = values_j[0] - values_i[0]
            dy = values_j[1] - values_i[1]
            plt.plot([x, x+dx], [y, y+dy], **kwargs)

            j += 1
        i += 1

    # label the cluster
    text = "$\\alpha^{" + str(orbit_index) + "}$"
    plt.annotate(text, label_xy, **cluster_annotate_kwargs)


# Get cluster sites
clust = None
with open(clust_json_path, 'r') as f:
    clust = json.load(f)


# Plot cluster sites
cluster_colorcycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# x-y
v1 = np.array([1.0, 0.0, 0.0]) # horizontal
v2 = np.array([0.0, 1.0, 0.0]) # vertical
cabinet = None
view_basis = make_cartesian_view_basis(v1, v2)

# Create figure

eps = 1e-5
orbit_index = 0
curr_x = 1
curr_y = 0
curr_z = 0
max_x = 1
min_y = curr_y
row_max_clust_y = 0
curr_n_sites = 1
for orbit in clust['orbits']:

    sites = orbit['prototype']['sites']

    if len(sites) < 1:
        orbit_index += 1
        continue

    # Start a new row of clusters
    if len(sites) != curr_n_sites or curr_x >= x_stop:
        curr_x = 1
        curr_y = curr_y - row_max_clust_y
        row_max_clust_y = 0
        curr_n_sites = len(sites)

    # Pick cluster color
    color_i = orbit_index % len(cluster_colorcycle)
    cluster_color = cluster_colorcycle[color_i]

    # Make supercell containing the cluster
    T_clust, ijk_min = get_cluster_box_2d(sites)
    print(T_clust)
    print(ijk_min)

    total_shift_ijk = np.array([curr_x, curr_y, 0]) \
        - ijk_min + np.array([0, -T_clust[1][1], 0], dtype=int)
    plot_cluster(view_basis, L, prim_cart_basis, orbit_index, sites,
        cluster_shift_ijk=total_shift_ijk, cabinet=cabinet)

    curr_x += T_clust[0][0]
    if T_clust[1][1] > row_max_clust_y:
        row_max_clust_y = T_clust[1][1]

    if curr_x > max_x:
        max_x = curr_x
    if curr_y - row_max_clust_y < min_y:
        min_y = curr_y - row_max_clust_y

    orbit_index += 1
    print("---\n")

print("max_x:", max_x)
print("curr_y:", curr_y)

# Make room for unit cell
min_y -= 3

# Supercell to use
T_super = np.array([
    [max_x, 0, 0],
    [0, -min_y, 0],
    [0, 0, 1]], dtype=int)
background_shift_ijk = np.array([0, min_y, 0], dtype=int)
df = supercell_sites_df(prim, T_super, shift=background_shift_ijk, pad=True)
df = df[df['z'] < 0.95]
_df = view_sites(view_basis, df, asym_scatter_kwargs=asym_scatter_kwargs,
    cabinet=cabinet)
plot_background(view_basis, L, _df, T_super, background_shift_ijk,
    cabinet=cabinet, axis_labels=['a', 'a', 'c'], hex=True)


# # Show it!
# plt.show()

# # Save it!
fig = plt.gcf()
fig.set_size_inches(figure_inches_wide, figure_inches_tall)
fig.tight_layout()
plt.savefig("clusters_2d.pdf")
