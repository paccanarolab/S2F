from scipy import sparse
from scipy.sparse import linalg

from diffusion import Diffusion
from Utils import Utilities


class LabelWeightedPropagation(Diffusion):

    def __init__(self, graph, proteins, terms):
        r"""
        Label Weighted Propagation Algorithm

        Parameters
        ----------
        graph : scipy.sparse
            Graph in which the labels should be diffused
            (before the kernel is built)
        proteins : pandas.DataFrame
            Indices of the proteins that conform the graph.
            This DataFrame can be built using the
            stand-alone 'utils' command
        terms : pandas.DataFrame
            Indices of the GO terms that will be mapped to
            the diffused seed. This DataFrame can be built using the
            stand-alone 'utils' command
        """
        super(LabelWeightedPropagation, self).__init__()
        self.graph = graph
        self.proteins = proteins
        self.terms = terms
        self.kernel_params = {}
        self.latest_diffusion = None
        self.kernel = None
        self.initial_guess = None

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
        Diffusion._write_results(self.latest_diffusion,
                                 self.proteins,
                                 self.terms, filename)

    def diffuse(self, initial_guess, **kwargs):
        r"""
        Diffuses the initial labelling `initial_guess`
        into the built kernel, if the kernel hasn't been built, it will
        build it using `compute_kernel`
        Parameters
        ----------
        initial_guess : scipy.sparse matrix
            The initial labelling matrix that will be diffused on the
            graph, shapes must be consistent to the given graph.

        Returns
        -------
        scipy.sparse.coo_matrix
            the new labelling after performing the label propagation

        Notes
        -----
        This propagation method requires the kernel to b
        recalculated every time
        """
        self.tell('Got the initial labelling')
        self.initial_guess = initial_guess.tocsc()
        self.compute_kernel(**kwargs)
        self.tell('Starting diffusion...')

        self.latest_diffusion = self.kernel * self.initial_guess
        self.latest_diffusion = self.latest_diffusion.tocoo()
        self.tell('done')

        # iterative version, keeping just in case
        '''F = None
        n = self.initial_guess.shape[0]
        m = self.initial_guess.shape[1]
        for i in range(m):
            k = self.initial_guess[:, i].todense() +\
                self.kernel_params['gamma']
            k = k.reshape((1, n))
            K = sparse.spdiags(k, 0, n, n)
            kernel = linalg.inv((K + self.kernel).tocsc())*K
            if F is None:
                F = kernel * self.initial_guess[:, i]
            else:
                F = sparse.hstack([F, kernel * self.initial_guess[:, i]])
            if i % 10 == 0:
                self.tell('iteration {i}/{m}'.format(i=i, m=m))
        self.latest_diffusion = F.tocoo()
        self.tell('done')'''

        return self.latest_diffusion

    def compute_kernel(self, **kwargs):
        r"""
        Combines the S2F kernel with the initial labelling as weights.
        Check `S2FLabelPropagation` for a description on the S2F kernel.
        The description here takes once the
        Laplacian is already computed

        .. math:: (I + \lambda L)^{-1}

        where :math:`L` is the Laplacian of matrix :math:`W_{S2F}`

        Parameters
        ----------
        kwargs
            Parameters to compute the kernel, the following entries
            will be handled:
            * 'lambda' : float
                To be used in the transformation described above.
                Defaults to 1.0

        Notes
        -----
        The kernel in this diffusion method is related to the initial
        labelling, and therefore this method will be called from the
        `diffuse` method once, and the final kernel will be built
        iterating over the columns.
        """
        if self.initial_guess is not None and self.set_kernel_params(**kwargs):
            # check that we have all we need
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
            lambdaL = L_S2F.multiply(self.kernel_params['lambda'])

            # compute the max of every protein
            k = self.initial_guess.max(axis=1).todense() +\
                self.kernel_params['gamma']
            k = k.reshape((1, n))
            K = sparse.spdiags(k, 0, n, n)
            self.kernel = linalg.inv((K + lambdaL).tocsc()) * K

        else:
            if self.initial_guess is None:
                self.tell('Waiting for the initial labelling to be' +
                          'setup before computing the kernel...')
            else:
                self.warning('Wrong parameters in Compute Kernel')
                raise ValueError('Wrong parameters in Compute Kernel')

    def set_kernel_params(self, **kwargs):
        self.tell('reading kernel parameters')
        self.kernel_params['lambda'] = kwargs.get('lambda', 1.0)
        self.kernel_params['gamma'] = kwargs.get('gamma', 0.1)
        self.tell(self.kernel_params)
        return True
