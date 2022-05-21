import math
import numpy as np
import scipy.linalg

class StrainConverter(object):
    """
    Converts between unrolled strain metric values (6-element or less vector
    representing a symmetric strain metric) and the deformation tensor, F.

    L: lattice vector (column) matrix (3x3)
    F: deformation tensor (3x3)
    R: rotation matrix (3x3)
    U: right stretch tensor (3x3)
    I: identity matrix (3x3)
    E: symmetric strain metric (3x3)
    unrolled_E (standard basis): ['E_{xx}', 'E_{yy}', 'E_{zz}', '\sqrt(2)E_{yz}', '\sqrt(2)E_{xz}', '\sqrt(2)E_{xy}']

    L_strained = F @ L_ideal
    F = R @ U

    - Green-Lagrange strain metric ('GLstrain'): E = (1/2)*(F.transpose() @ F - I)
    - Hencky strain metric ('Hstrain'): E = (1/2)*ln(F.transpose() @ F)
    - Euler-Almansi strain metric ('EAstrain'): E = (1/2)*(I−(F @ F.transpose()).inverse())
    - Biot strain metric ('Bstrain'): E = U - I
    - Right stretch tensor ('Ustrain'): E = U


    Attributes
    ----------
    metric: str (optional, default='Ustrain')
            Choice of strain metric, one of: 'Ustrain', 'GLstrain', 'Hstrain', 'EAstrain', 'Bstrain'

    basis: array of shape (6, dim)
        Choice of basis for unrolled_E, in terms of the standard basis.

            unrolled_E_in_standard_basis = basis @ unrolled_E_in_this_basis

    dim: integer
        Dimension (number of columns) of basis

    """
    @staticmethod
    def symmetry_adapted_basis():
        """
        Returns symmetry adapted basis, B, common to all 3d point groups:
            B[:,0] = e_1 = [1/sqrt(3), 1/sqrt(3), 1/sqrt(3), 0, 0, 0]
            B[:,1] = e_2 = [1/sqrt(2), -1/sqrt(2), 0.0, 0, 0, 0]
            B[:,2] = e_3 = [-1/sqrt(6), -1/sqrt(6), 2/sqrt(6), 0, 0, 0]
            B[:,3] = e_4 = [0, 0, 0, 1, 0, 0]
            B[:,4] = e_5 = [0, 0, 0, 0, 1, 0]
            B[:,5] = e_6 = [0, 0, 0, 0, 0, 1]

        """
        return np.array([
            [1./math.sqrt(3.), 1./math.sqrt(3.), 1./math.sqrt(3.), 0.0, 0.0, 0.0],
            [1./math.sqrt(2.), -1./math.sqrt(2.), 0., 0.0, 0.0, 0.0],
            [-1./math.sqrt(6.), -1./math.sqrt(6.), 2./math.sqrt(6.), 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]).transpose()

    def __init__(self, metric='Ustrain', basis=None):
        """
        Arguments
        ---------
        metric: str (optional, default='Ustrain')
            Choice of strain metric, one of: 'Ustrain', 'GLstrain', 'Hstrain', 'EAstrain', 'Bstrain'

        basis: array-like of shape (6, dim)
            User-choice of basis for unrolled_E, in terms of the standard basis.

                unrolled_E_in_standard_basis = basis @ unrolled_E

        """

        if basis is None:
            basis = np.identity(6)
        if len(basis.shape) != 2:
            raise Exception("Error in StrainConverter: len(basis.shape) != 2")
        if basis.shape[0] != 6:
            raise Exception("Error in StrainConverter: basis.shape[0] != 6")
        if basis.shape[1] > 6:
            raise Exception("Error in StrainConverter: basis.shape[1] > 6")

        if metric == 'GLstrain':
            self.desc = 'Green-Lagrange strain metric'
        elif metric == 'Hstrain':
            self.desc = 'Hencky strain metric'
        elif metric == 'EAstrain':
            self.desc = 'Euler-Almansi strain metric'
        elif metric == 'Bstrain':
            self.desc = 'Biot strain metric'
        elif metric == 'Ustrain':
            self.desc = 'Right stretch tensor'
        else:
            raise Exception("Error in StrainConverter: Unrecognized strain metric: " + str(metric))

        self.metric = metric
        self.basis = basis
        self.dim = basis.shape[1]

        self._basis_pinv = np.linalg.pinv(basis)

    def _to_unrolled_E(self, rolled_E):
        """Convert 3x3 symettric matrix strain metric to dim-element unrolled vector value"""
        E = rolled_E
        w = math.sqrt(2.)
        unrolled_E_in_standard_basis = np.array([E[0,0], E[1,1], E[2,2], w*E[1,2], w*E[0,2], w*E[0,1]])
        return self._basis_pinv @ unrolled_E_in_standard_basis

    def _to_rolled_E(self, unrolled_E):
        """Convert dim-element unrolled vector strain metric to symmetric 3x3 matrix value"""
        e = np.matmul(self.basis, unrolled_E)
        w = math.sqrt(2.)
        return np.array([[   e[0], e[5]/w, e[4]/w ],
                         [ e[5]/w,   e[1], e[3]/w ],
                         [ e[4]/w, e[3]/w,   e[2] ]])

    def to_F(self, unrolled_E):
        """Convert from unrolled_E to F

        Arguments
        -------
        unrolled_E: 1d array of size self.dim
            Unrolled strain metric value in terms of self.basis.

        Returns
        ---------
        F: array of shape (3,3)
            The deformation tensor

        """
        E = self._to_rolled_E(unrolled_E)
        if self.metric == 'GLstrain':
            # E = (1/2)*(F.transpose()*F - I)
            return scipy.linalg.sqrtm(2.*E + np.identity(3))
        elif self.metric == 'Hstrain':
            # E = (1/2)*ln(F.transpose()*F)
            D, V = np.linalg.eigh(E)
            return V @ np.diag(np.exp(D)) @ V.transpose()
        elif self.metric == 'EAstrain':
            # E = (1/2)*(I−(F*F.transpose()).inverse())
            return np.linalg.inv(scipy.linalg.sqrtm(np.identity(3) - 2.*E))
        elif self.metric == 'Bstrain':
            # E = U - I
            return E + np.identity(3)
        elif self.metric == 'Ustrain':
            # E = U
            return E
        else:
            raise Exception("Unrecognized strain metric: " + str(self.metric))

    def from_F(self, F):
        """Convert from F to unrolled_E

        Arguments
        ---------
        F: array of shape (3,3)
            The deformation tensor

        Returns
        -------
        unrolled_E: 1d array of size self.dim
            Unrolled strain metric value in terms of self.basis.
        """
        C = F.transpose() @ F
        if self.metric == 'GLstrain':
            # E = (1/2)*(F.transpose()*F - I)
            E = (1./2.)*(C - np.identity(3))
        elif self.metric == 'Hstrain':
            # E = (1/2)*ln(F.transpose()*F)
            E = (1./2.)*scipy.linalg.logm(C)
        elif self.metric == 'EAstrain':
            # E = (1/2)*(I−(F*F.transpose()).inverse())
            E = (1./2.)*(np.identity(3) - np.linalg.inv(C))
        elif self.metric == 'Bstrain':
            # E = U - I
            E =  scipy.linalg.sqrtm(C) - np.identity(3)
        elif self.metric == 'Ustrain':
            # E = U
            E = scipy.linalg.sqrtm(C)
        else:
            raise Exception("Unrecognized strain metric: " + str(self.metric))
        return self._to_unrolled_E(E)
