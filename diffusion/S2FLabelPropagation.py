from Utils import *
from diffusion import Diffusion

from scipy import sparse
from scipy.sparse import linalg

import numpy as np

class S2FLabelPropagation(Diffusion):

    def __init__(self, graph, proteins, terms):
        r"""
        S2F Label Propagation Algorithm

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
        super(S2FLabelPropagation, self).__init__()
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
        self.latest_diffusion = self.kernel * initial_guess.todense()
        self.latest_diffusion = sparse.coo_matrix(self.latest_diffusion)
        self.tell('done')
        return self.latest_diffusion

    def compute_kernel(self, **kwargs):
        r"""
        Computes the kernel required by S2F making the following transformations:

        .. math:: J_{ij} = \frac{\sum_{k}W_{ik}W_{jk}}{\sum_{k}W_{ik} + \sum_{k}W_{jk} - \sum_{k}W_{ik}W_{jk}}

        .. math:: {W_{\text{S2F}}}_{ij} = \frac{1}{2}\left(\frac{1}{d_i}+ \frac{1}{d_j}\right)J_{ij}W_{ij}

        where :math:`W` is the graph that was provided to the instance. The final kernel is

        .. math:: (I + \lambda L)^{-1}

        where :math:`L` is the Laplacian of matrix :math:`W_{S2F}`

        Parameters
        ----------
        kwargs
            Parameters to compute the kernel, the following entries will be handled:
            * 'lambda' : float
                To be used in the transformation described above. Defaults to 1.0

        Notes
        -----
        The kernel is available in `self.kernel`
        """
        if self.set_kernel_params(**kwargs):
            self.tell('Diffusion Kernel computation started...')
            n = self.graph.shape[0]
            self.tell('Computing Jaccard matrix...')
            J = sparse.csc_matrix(Utilities.jaccard(self.graph.todense()))
            self.tell('Jaccard * graph...')
            JW = J.multiply(self.graph)

            self.tell('Normalisation...')
            degree = 1.0/JW.sum(axis=1)
            W_S2F = 0.5 * (JW.multiply(degree) + JW.multiply(degree.T))

            self.tell('Computing Laplacian...')
            degree = W_S2F.sum(axis=1)
            D_S2F = sparse.spdiags(degree.T, 0, n, n)
            L_S2F = D_S2F - W_S2F
            IlambdaL = sparse.identity(n) + L_S2F.multiply(self.kernel_params['lambda'])

            self.tell(r'Inverting (I + \lambda W_S2F)...')
            self.kernel = np.linalg.inv(IlambdaL.todense())
            self.tell('Kernel built')
        else:
            self.warning('Wrong parameters in Compute Kernel')
            raise ValueError('Wrong parameters in Compute Kernel')

    def set_kernel_params(self, **kwargs):
        self.tell('reading kernel parameters...')
        self.kernel_params['lambda'] = kwargs.get('lambda', 1.0)
        self.tell(self.kernel_params)
        return True
