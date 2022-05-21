import math
import numpy as np

def make_cartesian_view_basis(v1=[1., 0., 0.], v2=[0., 0., 1.]):
    """
        v1: array of size (3,)
            Horizontal vector in the plane.
        v2: array of size (3,)
            Vertical vector in the plane.

    """
    v1 = np.array(v1)
    v2 = np.array(v2)
    assert v1.shape[0] == 3
    assert v2.shape[0] == 3

    b1 = v1 / np.linalg.norm(v1)
    b2 = v2 - (v2 @ b1) * b1
    b2 = b2 / np.linalg.norm(b2)
    b3 = np.cross(b1, b2)

    return np.array([b1, b2, b3]).transpose()


def make_view_basis(lattice_column_matrix, indices_1=[1., 0., 0.],
                    indices_2=[0., 0., 1.]):
    """
        lattice_column_matrix: array of shape (3,3)
            Ideal lattice vectors, as a column vector matrix
        indices_1: array of size (3,)
            Miller indices for first vector between lattice points in the plane.
        indices_2: array of size (3,)
            Miller indices for second vector between lattice points in the
            plane.

    """
    L = np.array(lattice_column_matrix)
    i1 = np.array(indices_1)
    i2 = np.array(indices_2)
    v1 = L @ i1
    v2 = L @ i2

    v1 = np.array(v1)
    v2 = np.array(v2)
    assert v1.shape[0] == 3
    assert v2.shape[0] == 3

    b1 = v1 / np.linalg.norm(v1)
    b2 = v2 - (v2 @ b1) * b1
    b2 = b2 / np.linalg.norm(b2)
    b3 = np.cross(b1, b2)

    return np.array([b1, b2, b3]).transpose()
