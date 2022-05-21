import math
import matplotlib.pyplot as plt
import numpy as np
import casm_figures.view

def begin_frac_2d():
    return np.array([
        [0., 0., 0.],
        [0., 1., 0.],
        [0., 0., 0.],
        [1., 0., 0.]
    ]).transpose()

def end_frac_2d():
    return np.array([
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [1., 1., 0.]
    ]).transpose()

def begin_frac():
    return np.array([
        [0., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [0., 1., 1.],
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 0., 1.],
        [1., 0., 1.],
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [1., 1., 0.]
    ]).transpose()

def end_frac():
    return np.array([
        [1., 0., 0.],
        [1., 1., 0.],
        [1., 0., 1.],
        [1., 1., 1.],
        [0., 1., 0.],
        [1., 1., 0.],
        [0., 1., 1.],
        [1., 1., 1.],
        [0., 0., 1.],
        [1., 0., 1.],
        [0., 1., 1.],
        [1., 1., 1.]
    ]).transpose()

def hex_begin_frac():
    return np.array([
        [0., 0., 0.],  # center vertical
        [1., 0., 0.],  # bottom hex
        [1., 1., 0.],
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.],
        [1., 0., 0.],  # vertical-bottom
        [1., 1., 0.],
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.],
        [1., 0., 1.],   # top hex
        [1., 1., 1.],
        [0., 1., 1.],
        [-1., 0., 1.],
        [-1., -1., 1.],
        [0., -1., 1.]
    ]).transpose()

def hex_end_frac():
    return np.array([
        [0., 0., 1.],  # center vertical
        [1., 1., 0.],  # bottom hex
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.],
        [1., 0., 0.],
        [1., 0., 1.],  # vertical-top
        [1., 1., 1.],
        [0., 1., 1.],
        [-1., 0., 1.],
        [-1., -1., 1.],
        [0., -1., 1.],
        [1., 1., 1.],   # top hex
        [0., 1., 1.],
        [-1., 0., 1.],
        [-1., -1., 1.],
        [0., -1., 1.],
        [1., 0., 1.]
    ]).transpose()

def hex_begin_frac_2d():
    return np.array([
        [1., 0., 0.],  # bottom hex
        [1., 1., 0.],
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.]
    ]).transpose()

def hex_end_frac_2d():
    return np.array([
        [1., 1., 0.],  # bottom hex
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.],
        [1., 0., 0.]
    ]).transpose()

def corners_frac():
    return np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [1., 1., 0.],
        [0., 0., 1.],
        [1., 0., 1.],
        [0., 1., 1.],
        [1., 1., 1.]
    ]).transpose()

def corners_frac_2d():
    return np.array([
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [1., 1., 0.]
    ]).transpose()

def hex_corners_frac():
    return np.array([
        [1., 1., 0.],  # bottom hex
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.],
        [1., 0., 0.],
        [1., 1., 1.],   # top hex
        [0., 1., 1.],
        [-1., 0., 1.],
        [-1., -1., 1.],
        [0., -1., 1.],
        [1., 0., 1.]
    ]).transpose()

def hex_corners_frac_2d():
    return np.array([
        [1., 1., 0.],  # bottom hex
        [0., 1., 0.],
        [-1., 0., 0.],
        [-1., -1., 0.],
        [0., -1., 0.],
        [1., 0., 0.]
    ]).transpose()

def apply_cabinet(cabinet, values):
    """Apply 'cabinet' perspective to values. This is similar to how most people draw 3d shapes (like a kitchen cabinet) by hand.

        cabinet: (f, theta)
            f: factor indicating fraction of "real" length displayed for vectors perpendicular to the viewing plane
            theta: angle the vectors are displayed at
        values: array of shape (3, n)
    """
    if cabinet is not None:
        for l in range(values.shape[1]):
            values[0,l] += -cabinet[0] * values[2,l] * math.cos(cabinet[1])
            values[1,l] += -cabinet[0] * values[2,l] * math.sin(cabinet[1])

def lattice_points_2d(lattice_column_matrix, indices_1=[1., 0., 0.], indices_2=[0., 1., 0.], lrange=None):
    """
        lattice_column_matrix: array of shape (3,3)
            Ideal lattice vectors, as a column vector matrix
        indices_1: array of size (3,)
            Miller indices for first vector between lattice points in the plane.
        indices_2: array of size (3,)
            Miller indices for second vector between lattice points in the
            plane.
        lrange: The lattice points to include, formatted as
            `[axis_1_begin, axis_1_end, axis_2_begin, axis_2_end]`
    """
    L = np.array(lattice_column_matrix)
    i1 = np.array(indices_1)
    i2 = np.array(indices_2)
    v1 = L @ i1
    v2 = L @ i2
    p1 = np.linspace(lrange[0], lrange[1], lrange[1]-lrange[0]+1)
    p2 = np.linspace(lrange[2], lrange[3], lrange[3]-lrange[2]+1)
    pts = np.zeros((3, len(p1)*len(p2)))
    l = 0
    for i in range(len(p1)):
        for j in range(len(p2)):
            pts[:,l] = v1*p1[i] + v2*p2[j]
            l += 1
    return pts

def lattice_points_3d(lattice_column_matrix, indices_1=[1., 0., 0.], indices_2=[0., 1., 0.], indices_3=[0., 0., 1.], lrange=None):
    """
        lattice_column_matrix: array of shape (3,3)
            Ideal lattice vectors, as a column vector matrix
        indices_1: array of size (3,)
            Miller indices for first vector between lattice points in the plane.
        indices_2: array of size (3,)
            Miller indices for second vector between lattice points in the
            plane.
        indices_3: array of size (3,)
            Miller indices for third vector between lattice points.
        lrange: The lattice points to include, formatted as
            `[begin_1, end_1, begin_2, end_2, begin_3, end_3]`
    """
    L = np.array(lattice_column_matrix)
    i1 = np.array(indices_1)
    i2 = np.array(indices_2)
    i3 = np.array(indices_3)
    v1 = L @ i1
    v2 = L @ i2
    v3 = L @ i3
    p1 = np.linspace(lrange[0], lrange[1], lrange[1]-lrange[0]+1)
    p2 = np.linspace(lrange[2], lrange[3], lrange[3]-lrange[2]+1)
    p3 = np.linspace(lrange[4], lrange[5], lrange[5]-lrange[4]+1)
    pts = np.zeros((3, len(p1)*len(p2)*len(p3)))
    l = 0
    for i in range(len(p1)):
        for j in range(len(p2)):
            for k in range(len(p3)):
                pts[:,l] = v1*p1[i] + v2*p2[j] + v3*p3[k]
                l += 1
    return pts

def plot_cartesian_axes(view_basis, scale=1., arrow_kwargs=None,
    annotate=True, annotate_kwargs=None, cabinet=None,
    shift=None, dim=3):
    """
    Arguments:
        scale: Length of vectors to draw
        shift: If provided, vectors begin at the Cartesian coordinates `shift`.
        cabinet: Draw with "cabinet" perspective. Expect tuple (length_factor,
            angle), where length_factor is length of out-of-plane lines, and
            angle is the angle the lines are written at. Typical values would
            be (0.5, math.pi/6.).

    """
    if annotate_kwargs is None:
        annotate_kwargs = {
            'textcoords': 'offset points',
            'xytext': (3, 3)
            }
    plot_lattice_cell(view_basis, np.eye(3)*scale, center=False, shift=shift,
        cell=False, arrows=True, arrow_kwargs=arrow_kwargs, annotate=annotate,
        annotate_kwargs=annotate_kwargs, cabinet=cabinet,
        axis_labels=['x','y','z'], dim=dim)

def plot_lattice_change(view_basis, L_final, L_init, center=False,
    arrows=True, arrow_kwargs=None, cabinet=None, plot_kwargs=None,
    exaggeration_factor=1.0, hex=False, dim=3, shift_reference=False):
    """Plot arrows showing change in lattice cell corners

    Arguments:
        view_basis: columns are plot viewing axes (horizontal, vertical, out-of-
            plane)
        L_final: final lattice
        L_init: initial lattice
        center: If True, center lattice cells at (0., 0., 0.). If False, begin
            lattice vectors at (0., 0., 0.).
        arrows: If True, draw arrows for lattice vectors.
        cabinet: Draw with "cabinet" perspective. Expect tuple (length_factor,
            angle), where length_factor is length of out-of-plane lines, and
            angle is the angle the lines are written at. Typical values would
            be (0.5, math.pi/6.).
        exaggeration_factor: Exagerate length of lines showing change in the
            lattice cell corners. A value of 1.0 is no exaggeration. A value of
            2.0 makes the arrows twice as long as the "real" change in position.
        shift_reference: Set to True to shift reference lattice. Example: Use
            if center=False to get symmetric arrows for pure strain.

    """
    basis_inv = np.linalg.pinv(view_basis)

    _corners_frac = None
    if hex is True:
        if dim == 3:
            _corners_frac = hex_corners_frac()
        elif dim == 2:
            _corners_frac = hex_corners_frac_2d()
    else:
        if dim == 3:
            _corners_frac = corners_frac()
        elif dim ==2:
            _corners_frac = corners_frac_2d()

    cart_final = L_final @ _corners_frac
    cart_init = L_init @ _corners_frac

    values_init = basis_inv @ cart_init
    apply_cabinet(cabinet, values_init)
    values_final = basis_inv @ cart_final
    apply_cabinet(cabinet, values_final)

    # shift cell center to origin
    shift_frac = np.array([-0.5, -0.5, -0.5])
    if hex is True:
        shift_frac = np.array([0.0, 0.0, -0.5])
    for i in range(_corners_frac.shape[1]):
        _corners_frac[:,i] = _corners_frac[:,i] + shift_frac

    center_cart_final = L_final @ _corners_frac
    center_cart_init = L_init @ _corners_frac

    center_values_init = basis_inv @ center_cart_init
    apply_cabinet(cabinet, values_init)
    center_values_final = basis_inv @ center_cart_final
    apply_cabinet(cabinet, values_final)

    if plot_kwargs is None:
        plot_kwargs = {}
    f = exaggeration_factor
    for i in range(values_init.shape[1]):
        x = None
        y = None
        dx = None
        dy = None
        if center:
            x = center_values_init[0][i]
            y = center_values_init[1][i]
            dx = center_values_final[0][i] - x
            dy = center_values_final[1][i] - y
        elif shift_reference:
            x_final = values_final[0][i]
            y_final = values_final[1][i]
            dx = center_values_final[0][i] - center_values_init[0][i]
            dy = center_values_final[1][i] - center_values_init[1][i]
            x = x_final - dx
            y = y_final - dy

        else:
            x = values_init[0][i]
            y = values_init[1][i]
            dx = values_final[0][i] - x
            dy = values_final[1][i] - y
        # scaling: keep head in same place, lengthen arrow by f
        plt.arrow(x+dx-f*dx, y+dy-f*dy, f*dx, f*dy, **arrow_kwargs)


def plot_lattice_cell(view_basis, lattice_column_matrix, center=False,
    cell=True, arrows=True, arrow_kwargs=None, annotate=True,
    annotate_kwargs=None, cabinet=None, plot_kwargs=None,
    axis_labels=['a', 'b', 'c'], shift=None, hex=False, dim=3):
    """
    Arguments:
        cell: If True, draw all 12 cell edges
        plot_kwargs: Extra arguments for plt.plot, used drawing cell lines
        arrows: If True, draw arrows for lattice vectors.
        arrow_kwargs: Extra arguments for plt.arrow, used to draw lattice
            vectors
        annotate: If True, name lattice vectors
        axis_labels: Labels used for lattice vectors, if annotating
        annotate_kwargs: Extra arguments for plt.annotate, used to name lattice
            vectors
        cabinet: Draw with "cabinet" perspective. Expect tuple (length_factor,
            angle), where length_factor is length of out-of-plane lines, and
            angle is the angle the lines are written at. Typical values would
            be (0.5, math.pi/6.).
        center: If provided, draw lattice so that cell body center is located
            at the origin.
        shift: If provided, lattice vectors begin at the Cartesian
            coordinates `shift`. The argument `center` takes precedence over
            `shift`.
        hex: Draw using standard hexagonal cell
        dim: Use dim==3 to draw 3d lattice cells, and dim==2 to draw 2d lattice
            cells
    """
    basis_inv = np.linalg.pinv(view_basis)

    L = lattice_column_matrix
    L_inv = np.linalg.pinv(L)

    shift_frac = np.array([0.0, 0.0, 0.0])
    if center:
        shift_frac = np.array([-0.5, -0.5, -0.5])
        if hex is True:
            shift_frac = np.array([0.0, 0.0, -0.5])
    elif shift is not None:
        shift_frac = L_inv @ shift

    _begin_frac = None
    if hex is True:
        if dim == 3:
            _begin_frac = hex_begin_frac()
        elif dim == 2:
            _begin_frac = hex_begin_frac_2d()
    else:
        if dim == 3:
            _begin_frac = begin_frac()
        elif dim == 2:
            _begin_frac = begin_frac_2d()

    for i in range(_begin_frac.shape[1]):
        _begin_frac[:,i] = _begin_frac[:,i] + shift_frac

    _end_frac = None
    if hex is True:
        if dim == 3:
            _end_frac = hex_end_frac()
        elif dim == 2:
            _end_frac = hex_end_frac_2d()
    else:
        if dim == 3:
            _end_frac = end_frac()
        elif dim == 2:
            _end_frac = end_frac_2d()

    for i in range(_end_frac.shape[1]):
        _end_frac[:,i] = _end_frac[:,i] + shift_frac


    begin_cart = lattice_column_matrix @ _begin_frac
    end_cart = lattice_column_matrix @ _end_frac

    begin_values = basis_inv @ begin_cart
    apply_cabinet(cabinet, begin_values)
    end_values = basis_inv @ end_cart
    apply_cabinet(cabinet, end_values)

    if cell is True:
        if plot_kwargs is None:
            plot_kwargs = {}
        for i in range(begin_values.shape[1]):
            x = [begin_values[0][i], end_values[0][i]]
            y = [begin_values[1][i], end_values[1][i]]
            zorder = (begin_values[2][i] + end_values[2][i])/2.
            plt.plot(x, y, zorder=zorder, **plot_kwargs)

    # drawing and annotating arrows
    begin = np.zeros((3,3))
    for i in range(3):
        begin[:,i] = L @ shift_frac
    begin_values = basis_inv @ begin
    apply_cabinet(cabinet, begin_values)
    end_values = basis_inv @ (begin + lattice_column_matrix)
    apply_cabinet(cabinet, end_values)
    x = begin_values[0,:]
    y = begin_values[1,:]
    z = begin_values[2,:]
    end_x = end_values[0,:]
    end_y = end_values[1,:]
    end_z = end_values[2,:]
    dx = end_values[0,:] - begin_values[0,:]
    dy = end_values[1,:] - begin_values[1,:]
    dz = end_values[2,:] - begin_values[2,:]

    if arrows is True:
        if arrow_kwargs is None:
            arrow_kwargs = {}
        for i in range(dim):
            if math.isclose(dx[i],0.0) and math.isclose(dy[i], 0.0):
                # marker=r'$\odot$'
                # if dz[i] < 0.0:
                #     marker=r'$\bigotimes$'
                # # plt.scatter(x[i], y[i], s=150, marker='o', facecolor='white')
                # plt.scatter(x[i], y[i], s=150, marker=marker, facecolor='white',
                #     edgecolors=arrow_kwargs.get('color','black'))
                continue
            else:
                zorder = (z[i] + (z[i] + dz[i]))/2.
                plt.arrow(x[i], y[i], dx[i], dy[i], zorder=zorder+0.1,
                    **arrow_kwargs)

    if annotate is True:
        if annotate_kwargs is None:
            annotate_kwargs = {
                'textcoords': 'offset points',
                'xytext': (3, 3)
                }
        for i in range(dim):
            if math.isclose(dx[i],0.0) and math.isclose(dy[i], 0.0):
                continue
            plt.annotate(axis_labels[i], (end_x[i], end_y[i]),
                         zorder=end_z[i]+0.1, **annotate_kwargs)

def plot_lattice_cell_2d(view_basis, lattice_column_matrix, arrows=True,
    arrow_kwargs=None, annotate=True, annotate_kwargs=None, plot_kwargs=None):
    basis_inv = np.linalg.pinv(view_basis)
    begin_cart = lattice_column_matrix @ begin_frac_2d()
    end_cart = lattice_column_matrix @ end_frac_2d()

    begin_values = basis_inv @ begin_cart
    end_values = basis_inv @ end_cart

    if plot_kwargs is None:
        plot_kwargs = {}
    for i in range(begin_values.shape[1]):
        x = [begin_values[0][i], end_values[0][i]]
        y = [begin_values[1][i], end_values[1][i]]
        plt.plot(x, y, **plot_kwargs)

    values = basis_inv @ lattice_column_matrix
    if arrows is True:
        if arrow_kwargs is None:
            arrow_kwargs = {}
        plt.arrow(0, 0, values[0][0], values[1][0], **arrow_kwargs)
        plt.arrow(0, 0, values[0][1], values[1][1], **arrow_kwargs)

    if annotate is True:
        if annotate_kwargs is None:
            annotate_kwargs = {
                'textcoords': 'offset points',
                'xytext': (3, 3)
                }
        plt.annotate("a", (values[0][0], values[1][0]), **annotate_kwargs)
        plt.annotate("b", (values[0][1], values[1][1]), **annotate_kwargs)

def plot_points(view_basis, points, cabinet=None, **kwargs):
    basis_inv = np.linalg.pinv(view_basis)
    values = basis_inv @ points
    apply_cabinet(cabinet, values)
    plt.scatter(values[0,:], values[1,:], **kwargs)
