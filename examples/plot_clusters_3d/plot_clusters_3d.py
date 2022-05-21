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
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'pdf.fonttype': 42})
plt.rcParams.update({'ps.fonttype': 42})

root = "fcc_binary"
clust_json_path = "fcc_binary/basis_sets/bset.default/clust.json"
output_dirpath = "clusters_3d"
output_format = "pdf"
figsize=(10,10)

# Which supercell to use, S = L @ T, where S and L are the supercell and prim lattice vectors as columns of a 3x3 matrix, and T is the integer transformation matrix to generate the supercell
#T_conventional = np.eye(3)
T_conventional = np.array([
    [-1., 1., 1.],
    [1., -1., 1.],
    [1., 1., -1.]])

dim = 3

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

# Show conventional hexagonal cell
hex = False

# Set plotting parameters for sites by asymmetric unit
s_scale = 100
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
    cell_cell = True
    cell_plot_kwargs = {
        'color': 'black',
        'linestyle': 'dashed'}
    cell_arrows = True
    cell_arrow_kwargs = {
        'head_width':.05,
        'length_includes_head':True,
        'color':'black'
    }
    cell_annotate = True
    cell_annotate_kwargs = {
        'textcoords': 'offset points',
        'xytext': (-20, -20)
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
                      cell=cell_cell,
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
                      cell=cell_cell,
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
    cartesian_axes_scale = 0.5
    cartesian_axes_head_width = cartesian_axes_scale / 10.
    cartesian_axes_scale_origin = L @ (supercell_shift_ijk - np.array([0.5, 0.5, 0.5]))
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
    cluster_scatter_kwargs=None, cluster_line_kwargs=None,
    default_cluster_options=None, cluster_options={},
    cluster_annotate_kwargs=None, zorder=-1000, cabinet=None):

    view_basis_inv = np.linalg.pinv(view_basis)
    eps = 1e-5
    if cluster_annotate_kwargs is None:
        cluster_annotate_kwargs = {
            'textcoords': 'offset points',
            'xytext': (-25, 5),
            'zorder': 1000
            }
    if default_cluster_options is None:
        default_cluster_options = {
            'element': 0,
            'shift_ijk': [0, 0, 0]
        }
    if cluster_line_kwargs is None:
        cluster_line_kwargs = {
            'linewidth': 6
        }
    if cluster_scatter_kwargs is None:
        cluster_scatter_kwargs = {
            's': 300
        }
    label_xy = None

    i = 0
    while i < len(sites):
        site_i = sites[i]
        ijk_i = np.array(site_i[-3:], dtype=int) #+ cluster_shift_ijk
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
            ijk_j = np.array(site_j[-3:], dtype=int) #+ cluster_shift_ijk
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

# Create directory for writing figures
dirpath = copy.copy(output_dirpath)
i = 0
while os.path.exists(dirpath):
    i += 1
    dirpath = output_dirpath + "." + str(i)
output_dirpath = copy.copy(dirpath)
os.mkdir(output_dirpath)

# Begin drawing figures, one for each orbit
eps = 1e-5
orbit_index = 0
orbit_included_index = 0
for orbit in clust['orbits']:

    sites = orbit['prototype']['sites']

    if len(sites) < 1:
        orbit_index += 1
        continue

    # Pick cluster color
    color_i = orbit_included_index % len(cluster_colorcycle)
    cluster_color = cluster_colorcycle[color_i]

    # Plot sites, cell, cartesian axes (technically it's the foreground)
    plt.subplots(2,2, figsize=figsize)

    # Make supercell containing the cluster
    T, supercell_shift_ijk = get_cluster_box(L, prim_cart_basis, sites,
        T_conventional=T_conventional)
    print(sites)
    print(T, supercell_shift_ijk)
    print("----\n")
    df = supercell_sites_df(prim, T, shift=supercell_shift_ijk, pad=True)

    # 'cabinet'
    plt.subplot(2,2,2)
    v1 = np.array([1.0, 0.0, 0.0]) # horizontal
    v2 = np.array([0.0, 0.0, 1.0]) # vertical
    cabinet = (0.2, math.pi/6.)
    view_basis = make_cartesian_view_basis(v1, v2)
    _df = view_sites(view_basis, df, asym_scatter_kwargs=asym_scatter_kwargs,
        cabinet=cabinet)
    plot_background(view_basis, L, _df, T, supercell_shift_ijk,
        cabinet=cabinet, T_conventional=T_conventional, hex=hex)
    plot_cluster(view_basis, L, prim_cart_basis, orbit_index, sites,
        cabinet=cabinet)

    # x-y
    plt.subplot(2,2,1)
    v1 = np.array([1.0, 0.0, 0.0]) # horizontal
    v2 = np.array([0.0, 1.0, 0.0]) # vertical
    cabinet = None
    view_basis = make_cartesian_view_basis(v1, v2)
    _df = view_sites(view_basis, df, asym_scatter_kwargs=asym_scatter_kwargs,
        cabinet=cabinet)
    plot_background(view_basis, L, _df, T, supercell_shift_ijk,
        cabinet=cabinet, T_conventional=T_conventional, hex=hex)
    plot_cluster(view_basis, L, prim_cart_basis, orbit_index, sites,
        cabinet=cabinet)

    # x-z
    plt.subplot(2,2,3)
    v1 = np.array([1.0, 0.0, 0.0]) # horizontal
    v2 = np.array([0.0, 0.0, 1.0]) # vertical
    cabinet = None
    view_basis = make_cartesian_view_basis(v1, v2)
    _df = view_sites(view_basis, df, asym_scatter_kwargs=asym_scatter_kwargs,
        cabinet=cabinet)
    plot_background(view_basis, L, _df, T, supercell_shift_ijk,
        cabinet=cabinet, T_conventional=T_conventional, hex=hex)
    plot_cluster(view_basis, L, prim_cart_basis, orbit_index, sites,
        cabinet=cabinet)

    # y-z
    plt.subplot(2,2,4)
    v1 = np.array([0.0, 1.0, 0.0]) # horizontal
    v2 = np.array([0.0, 0.0, 1.0]) # vertical
    cabinet = None
    view_basis = make_cartesian_view_basis(v1, v2)
    _df = view_sites(view_basis, df, asym_scatter_kwargs=asym_scatter_kwargs,
        cabinet=cabinet)
    plot_background(view_basis, L, _df, T, supercell_shift_ijk,
        cabinet=cabinet, T_conventional=T_conventional, hex=hex)
    plot_cluster(view_basis, L, prim_cart_basis, orbit_index, sites,
        cabinet=cabinet)

    # Show it!
    # plt.show()

    # Save it!
    filename = "orbit." + str(orbit_index) + "." + output_format
    plt.savefig(os.path.join(output_dirpath, filename))

    orbit_index += 1
    orbit_included_index += 1
