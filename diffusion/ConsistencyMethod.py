from Utils import *
from diffusion import Diffusion

from scipy import sparse
from scipy.sparse import linalg

import numpy as np
import pandas as pd

class ConsistencyMethod(Diffusion):

    def __init__(self, graph, proteins, terms):
        r"""
        Consistency Method Label Propagation Algorithm

        Parameters
        ----------
        graph : scipy.sparse
            Graph in which the labels should be diffused (before the kernel is built)
        proteins : pandas.DataFrame
            Indices of the proteins that conform the graph. This DataFrame can be built using the
            stand-alone 'utils' command
        terms : pandas.DataFrame
            Indices of the GO terms that will be mapped to the diffused seed. This DataFrame can be built using the
            stand-alone 'utils' command
        """
        super(ConsistencyMethod, self).__init__()
        self.graph = graph
        self.proteins = proteins
        self.terms = terms
        self.kernel_params = {}
        self.latest_diffusion = None
        self.kernel = None

    def write_results(self, filename):
        r"""
        Write the results of the diffusion to the path pointed by `filename`
        the format will be TSV with the following columns:
            * protein
            * goterm
            * score

        Parameters
        ----------
        filename : str
            Path to write the results
        """
        Diffusion._write_results(self.latest_diffusion, self.proteins, self.terms, filename)

    def diffuse(self, initial_guess, **kwargs):
        r"""
        Diffuses the initial labelling `initial_guess` into the built kernel, if the kernel hasn't been built, it will
        build it using `compute_kernel`
        Parameters
        ----------
        initial_guess : scipy.sparse matrix
            The initial labelling matrix that will be diffused on the graph, shapes must be consistent to the given
            graph.

        Returns
        -------
        scipy.sparse.coo_matrix
            the new labelling after performing the label propagation

        Notes
        -----
        The final labelling is kept in `self.latest_diffusion`, for access convenience. This enables a subsequent
        call to `write_results` that does not require a re-calculation of the final labelling.
        """
        self.tell('Starting diffusion...')
        self.latest_diffusion = self.kernel * initial_guess
        self.latest_diffusion = self.latest_diffusion.tocoo()
        self.tell('done')
        return self.latest_diffusion

    def compute_kernel(self, **kwargs):
        r"""

        .. math:: (I - \alpha S)^{-1}

        Parameters
        ----------
        kwargs
            Parameters to compute the kernel, the following entries will be handled:
            * 'alpha' : float
                To be used in the transformation described above. Defaults to 1.0

        Notes
        -----
        The kernel is available in `self.kernel`
        """

        if self.set_kernel_params(**kwargs):

            self.tell('Diffusion Kernel computation started...')
            # build D
            sums = self.graph.sum(1)  # sum every column per row

            # in case there are unconnected parts in the matrix
            indNonZeros = np.where(sums != 0)
            diagonalValues = np.zeros(sums.shape)

            diagonalValues[indNonZeros] = 1.0 / np.sqrt(sums[indNonZeros])
            D = sparse.spdiags(diagonalValues.T, 0, diagonalValues.shape[0], diagonalValues.shape[0])

            # build S
            S = D * self.graph * D

            IalphaS = np.eye(S.shape[0]) - self.kernel_params['alpha'] * S

            self.tell(r'Inverting (I - \alpha S)...')
            self.kernel = linalg.inv(IalphaS.tocsc())
            self.tell('Kernel built')
        else:
            self.warning('Wrong parameters in Compute Kernel')
            raise ValueError('Wrong parameters in Compute Kernel')

    def set_kernel_params(self, **kwargs):
        self.tell('reading kernel parameters...')
        self.kernel_params['alpha'] = kwargs.get('alpha', 1.0)
        self.tell(self.kernel_params)
        return True
