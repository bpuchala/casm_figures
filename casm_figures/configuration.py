import copy
import math
import numpy as np
import pandas as pd
from casm.info import get_prim_info, get_supercell_info

def configuration_sites_df(prim, config, shift=[0, 0, 0], pad=True,
    hex_lim=None):
    """Create a Dataframe containing sites within a supercell

    Arguments:
        prim: Prim, as dict
        config: Configuration, as dict. (i.e. config.json or query "config"
            output)
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
    T = np.array(config["transformation_matrix_to_supercell"], dtype=int)

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

    # add dofs
    N = len(info["integral_site_coordinates"])
    info["dofs"] = []
    local_dof_dim = {}
    print(config["dof"]["occ"])
    for i in range(N):
        info["dofs"].append({'occ': config["dof"]["occ"][i]})
    print(info["dofs"])
    if "local_dofs" in config["dof"]:
        for key in config["dof"]["local_dofs"]:
            values = config["dof"]["local_dofs"][key]["values"]
            for i in range(N):
                for j in range(len(values[i])):
                    info["dofs"][i][key + "_" + str(j)] = values[i][j]

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
                    info["dofs"].append(info["dofs"][i])

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
        new_dofs = []
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
                new_dofs.append(info["dofs"][i])
        info["integral_site_coordinates"] = new_bijk
        info["frac_coordinate"] = new_frac
        info["dofs"] = new_dofs

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

    data = {
        "x":supercell_cart_basis[0,:],
        "y":supercell_cart_basis[1,:],
        "z":supercell_cart_basis[2,:],
        "asym":asym,
        "b":integral_site_coordinates[0,:],
        "i":integral_site_coordinates[1,:],
        "j":integral_site_coordinates[2,:],
        "k":integral_site_coordinates[3,:]}

    # add DoF
    for key in info["dofs"][0]:
        values = []
        for i in range(N):
            values.append(info["dofs"][i][key])
        data[key] = values

    # add species_name and species_unique_name
    species_name = []
    species_unique_name = []
    for i in range(N):
        b = data["b"][i]
        occ = data["occ"][i]
        unique_name = prim["basis"][b]["occupants"][occ]
        name = unique_name
        if "species" in prim and unique_name in prim["species"]:
            name = prim["species"][unique_name]["name"]
        species_name.append(name)
        species_unique_name.append(unique_name)
    data["species_name"] = species_name
    data["species_unique_name"] = species_unique_name

    return pd.DataFrame(data)
