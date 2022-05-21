import math
import numpy as np

def cluster_cart_coordinates(prim_L, prim_cart_basis, cluster_integral_sites):
    """Return cluster sites as columns of an array, in Cartesian coordinates
    """
    C = np.zeros((3, len(cluster_integral_sites)))
    i = 0
    for site in cluster_integral_sites:
        C[:,i] = prim_cart_basis[site[0]] + prim_L @ np.array(site[-3:])
        i += 1
    return C

def select_cluster_sites(supercell_sites_df, cluster_integral_sites):
    selected = []
    count = 0
    for row in supercell_sites_df.itertuples():
        scel_site = [row.b, row.i, row.j, row.k]
        in_cluster = False
        for site in cluster_integral_sites:
            if site == scel_site:
                in_cluster = True
                count += 1
                break
        selected.append(in_cluster)
    supercell_sites_df['selected'] = selected
    if count != len(cluster_integral_sites):
        raise Exception("Error in select_cluster_sites: Not all sites found")
    return supercell_sites_df

def get_cluster_box_2d(sites):
    """
    Arguments:
        sites: list of [b, i, j, k], for each site in cluster

    Returns:
        (T, ijk_min):
            T: supercell transformation matrix that includes cluster sites
            ijk_min: coordinates of the lower-left corner of T
    """
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
    """ Find a small supercell that contains the cluster

    Arguments:
        L: prim lattice column matrix
        prim_cart_basis: array of shape (3, n_basis), with Cartesian coordinates
            of the prim basis sites
        sites: list of [b, i, j, k], for each site in cluster
        T_conventional: default = np.eye(3, dtype=int)
            Alternate "conventional" supercell that should be used

    Returns:
        (T, ijk_min):
            T: supercell transformation matrix that, when pad=True, includes
                cluster sites, when multiplied L @ T_conventional @ T
            ijk_min: prim lattice frac coordinates of the lower-left corner of
                box
    """
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
