import copy
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from casm.info import get_prim_info, get_supercell_info
from casm_figures.lattice import plot_lattice_cell, plot_cartesian_axes, \
    apply_cabinet

def supercell_sites_df(prim, T, shift=[0, 0, 0], pad=True, hex_lim=None):
    """Create a Dataframe containing sites within a supercell

    Arguments:
        prim: Prim, as dict
        T: 3x3 integer transformation_matrix_to_super
        shift: size=3 integer ijk coordinate shift to apply to site
            coordinates. This is the number of prim unit cells to shift the
            sites.
        pad: include site images on cell boundaries (i.e. if there is a
            site at one corner, include a site on all eight corners)
        hex_lim: If not None, remove sites outside the standard prim hexagonal
            cell boundaries, centered at the origin. Expected value is the
            tuple (n_a, n_c), indicating the size of the hexagonal cell in
            'a' and 'c' lengths included. When using this option,
            make sure that T and shift are sufficient to generate enough
            sites that the hexagonal cell is filled. Assumes lattice has
            vectors: L = [[a,0.,0.],[-a/2.,sqrt(3.)/2.,0.],[0.,0.,c]]
            For example, to generate the standard hex cell:
                T = [[2.,0.,0.],[0.,2.,0.],[0.,0.,0.]]
                shift = [-1,-1,0]
                pad = True
                hex_lim = (1, 1)

    Returns:
        DataFrame, 1 row per supercell site, and columns:
            x: site x Cartesian coordinate
            y: site y Cartesian coordinate
            z: site z Cartesian coordinate
            asym: asymmetric unit index (indicates orbit of equivalent sites)
            b: sublattice index
            i: integer unit cell index, corresponding to prim lattice vector a
            j: integer unit cell index, corresponding to prim lattice vector b
            k: integer unit cell index, corresponding to prim lattice vector c

    """
    T = np.array(T, dtype=int)

    # shift (in prim ijk)
    shift = np.array(shift, dtype=int)
    L = np.array(prim["lattice_vectors"]).transpose()
    L_inv = np.linalg.pinv(L)

    # supercell lattice vectors
    S = L @ T

    # shift in fractional coordinates (of supercell lattice vectors)
    S_inv = np.linalg.pinv(S)
    shift_frac = S_inv @ (L @ shift)

    info = get_prim_info(prim=prim, properties=["asymmetric_unit"])

    sublattice_to_asym = {}
    asym = 0
    for orbit in info["asymmetric_unit"]:
        for b in orbit:
            sublattice_to_asym[b] = asym
        asym += 1

    # get prim cart coordinates
    I = np.eye(3, dtype=int)
    info = get_supercell_info(prim=prim,
                              transformation_matrix_to_super=I.tolist(),
                              properties=["cart_coordinate"])
    prim_cart_basis = np.array(info["cart_coordinate"]).transpose()

    # get supercell integral site coordinates
    info = get_supercell_info(prim=prim,
                              transformation_matrix_to_super=T.tolist(),
                              properties=[
                                    "integral_site_coordinates",
                                    "frac_coordinate"])

    # CASM output of supercell sites are not all "within"
    # move sites within supercell
    within = True
    if within:
        N = len(info["integral_site_coordinates"])
        for j in [0, 1, 2]:
            for i in range(N):
                frac = info["frac_coordinate"][i]
                bijk = info["integral_site_coordinates"][i]
                while math.isclose(frac[j], 1.0) or frac[j] > 1.0:
                    frac[j] -= 1.0
                    bijk[1] -= T[0][j]
                    bijk[2] -= T[1][j]
                    bijk[3] -= T[2][j]

                while not math.isclose(frac[j], 0.0) and frac[j] < 0.0:
                    frac[j] += 1.0
                    bijk[1] += T[0][j]
                    bijk[2] += T[1][j]
                    bijk[3] += T[2][j]

    # add a site at frac coordinate 1.0 for all sites at frac coordiante 0.0
    if pad:
        for j in [0, 1, 2]:
            N = len(info["integral_site_coordinates"])
            for i in range(N):
                frac = copy.copy(info["frac_coordinate"][i])
                if math.isclose(frac[j], 0.):
                    new_bijk = copy.copy(info["integral_site_coordinates"][i])
                    new_bijk[1] += T[0][j]
                    new_bijk[2] += T[1][j]
                    new_bijk[3] += T[2][j]
                    new_frac = copy.copy(frac)
                    new_frac[j] = 1.0
                    info["integral_site_coordinates"].append(new_bijk)
                    info["frac_coordinate"].append(new_frac)

    # apply shift
    if True:
        N = len(info["integral_site_coordinates"])
        for i in range(N):
            for j in range(3):
                info["integral_site_coordinates"][i][j+1] += shift[j]
                info["frac_coordinate"][i][j] += shift_frac[j]

    # optionally, trim sites outside hexagonal cell
    if hex_lim is not None:
        new_bijk = []
        new_frac = []
        N = len(info["integral_site_coordinates"])
        eps = 1e-5
        for i in range(N):
            keep = True
            site = info["integral_site_coordinates"][i]
            b = site[0]
            ijk = site[-3:]
            cart = prim_cart_basis[:,b] + L @ ijk
            frac = L_inv @ cart # this is prim lattice frac
            if frac[0] > hex_lim[0] + eps:
                keep = False
            if frac[0] < -hex_lim[0] - eps:
                keep = False
            if frac[1] < -hex_lim[0] + frac[0] - eps:
                keep = False
            if frac[1] < -hex_lim[0] - eps:
                keep = False
            if frac[1] > hex_lim[0] + eps:
                keep = False
            if frac[1] > hex_lim[0] + frac[0] + eps:
                keep = False
            if frac[2] < -eps:
                keep = False
            if frac[2] > hex_lim[1] + eps:
                keep = False
            if keep:
                new_bijk.append(info["integral_site_coordinates"][i])
                new_frac.append(frac)
        info["integral_site_coordinates"] = new_bijk
        info["frac_coordinate"] = new_frac

    # generate data frame
    N = len(info["integral_site_coordinates"])
    integral_site_coordinates = np.array(info["integral_site_coordinates"], dtype=int).transpose()
    supercell_cart_basis = np.zeros((3, N))
    asym = []
    for i in range(N):
        site = integral_site_coordinates[:,i]
        b = site[0]
        ijk = site[-3:]
        asym.append(sublattice_to_asym[b])
        supercell_cart_basis[:,i] = prim_cart_basis[:,b] + L @ ijk

    return pd.DataFrame({
        "x":supercell_cart_basis[0,:],
        "y":supercell_cart_basis[1,:],
        "z":supercell_cart_basis[2,:],
        "asym":asym,
        "b":integral_site_coordinates[0,:],
        "i":integral_site_coordinates[1,:],
        "j":integral_site_coordinates[2,:],
        "k":integral_site_coordinates[3,:]})

def view_sites(view_basis, _df, key="asym", scatter_kwargs={}, cabinet=None):
    """Create Dataframe containing coordinates in viewing basis and columns containing plt.scatter properties

    Arguments:
        view_basis: Basis used for plotting
        _df: DataFrame containing sites as rows
        key: Name of column to be used to set scatter properties
        scatter_kwargs: Scatter plot properties
    """

    colorcycle = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    # for provided types, set a default color if it does not exist
    t_index = 0
    for t in scatter_kwargs:
        if 'c' not in scatter_kwargs[t]:
            color_i = t_index % len(colorcycle)
            scatter_kwargs[t]['c'] = colorcycle[color_i]
            t_index += 1

    # Example asym_scatter_kwargs:
    # Dict of asymmetric unit index (str) -> plt.scatter properties
    # asym_scatter_kwargs = {
    #     '0': {
    #         'c': 'red',
    #         's': 40.0,
    #         'edgecolors': 'black'
    #         'marker': '$\\oplus$'
    #     },
    #     '1': {
    #         'c': 'blue',
    #         's': 40.0,
    #         'edgecolors': 'black'
    #     },
    #     '2': {
    #         'c': 'white',
    #         's': 80.0,
    #         'edgecolors': 'black'
    #     }
    # }

    df = _df.copy()
    view_basis_inv = np.linalg.pinv(view_basis)
    values = view_basis_inv @ df[['x','y','z']].values.transpose()
    apply_cabinet(cabinet, values)
    df['px'] = values[0,:]
    df['py'] = values[1,:]
    df['pz'] = values[2,:]
    df = df.sort_values('pz')
    print([val for val in df[key]])
    df['c'] = [scatter_kwargs[str(val)]['c'] for val in df[key]]
    df['s'] = [scatter_kwargs[str(val)].get('s', 40.0) for val in df[key]]
    df['edgecolors'] = [scatter_kwargs[str(val)].get('edgecolors', 'black') for val in df[key]]
    df['marker'] = [scatter_kwargs[str(val)].get('marker', 'o') for val in df[key]]
    return df

def plot_background(view_basis, L, df, T,
    T_conventional=None, unitcell=True, unitcell_kwargs={},
    supercell=True, supercell_kwargs={}, cartesian_axes=True,
    cartesian_axes_kwargs={}):

    # Typical kwargs
    # cell_kwargs = {
    #     'plot_kwargs': {
    #         'color': 'black',
    #         'linestyle': 'dashed'
    #     },
    #     'arrow_kwargs': {
    #         'head_width':.1,
    #         'length_includes_head':True,
    #         'color':'black'
    #     },
    #     'annotate_kwargs': {
    #         'textcoords': 'offset points',
    #         'xytext': (-10, -10)
    #     }
    # }
    # cartesian_axes_kwargs = {
    #     'scale': 1.0,
    #     'shift': np.array([-2, 0, 0]),
    #     'arrow_kwargs': {
    #         'head_width': 0.1,
    #         'color': 'black',
    #         'head_width':cartesian_axes_head_width,
    #         'length_includes_head':True
    #     }
    # }

    if T_conventional is None:
        T_conventional = np.eye(3, dtype=int)

    # plot unit cell
    if unitcell:
        kwargs = copy.copy(unitcell_kwargs)
        plot_lattice_cell(view_basis, L @ T_conventional, **kwargs)
    # plot super cell
    if supercell:
        kwargs = copy.copy(supercell_kwargs)
        plot_lattice_cell(view_basis, L @ T, **kwargs)

    # Plot sites
    N = len(df['px'])
    for i in range(N):
        plt.scatter(x=[df['px'].values[i]], y=[df['py'].values[i]],
            c=[df['c'].values[i]], s=[df['s'].values[i]],
            edgecolors=df['edgecolors'].values[i], linewidths=1,
            zorder=df['pz'].values[i])

    # Plot Cartesian axes
    if cartesian_axes:
        kwargs = copy.copy(cartesian_axes_kwargs)
        plot_cartesian_axes(view_basis, **kwargs)
