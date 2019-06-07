from Utils import *
from graphs import Graph

from scipy import sparse
from scipy.sparse.linalg import lsqr
from sklearn.preprocessing import binarize

import pandas as pd
import numpy as np
import os


class Combination(Graph):

    def __init__(self, proteins, collection, homology, seed):
        super(Combination, self).__init__()
        self.combined = None
        self.proteins = proteins
        self.collection = collection
        self.homology = homology
        self.seed = seed

    def get_graph(self, **kwargs):
        return Graph.to_sparse_matrix(self.combined.T)

    def write_graph(self, filename):
        g = self.get_graph()
        data = {
            'protein1': g.row,
            'protein2': g.col,
            'score': g.data
        }
        combined_df = pd.DataFrame(data)
        combined_df = combined_df.merge(self.proteins.reset_index(), left_on='protein1', right_on='protein idx')
        combined_df = combined_df.merge(
            self.proteins.reset_index(), left_on='protein2',
            right_on='protein idx', suffixes=['1', '2'])
        combined_df[['protein id1', 'protein id2', 'score']].to_csv(filename, sep='\t', index=False, header=None)

    def compute_graph(self):
        self.tell('Building the kernel matrix...')
        kernels = None
        kernel_order = []
        for k, g in self.collection.items():
            self.tell('adding', k, 'graph...')
            kernel_order.append(k)
            if kernels is None:
                kernels = Graph.to_sparse_vector(sparse.triu(g, 1).tocoo())
            else:
                kernels = sparse.hstack([kernels, Graph.to_sparse_vector(sparse.triu(g, 1).tocoo())])
        self.tell('adding homology graph...')
        kernels = sparse.hstack([kernels, Graph.to_sparse_vector(sparse.triu(self.homology, 1).tocoo())]).tocsc()

        # normalise using a for loop because for some reason doing this with a matrix expression uses a different method
        # with more precision errors
        self.tell('normalising kernels')
        k_norm = kernels[:, 0].multiply(1 / (kernels[:, 0].max(axis=0).toarray() + np.finfo(float).eps))
        for i in range(1, kernels.shape[1]):
            inv_max = 1 / (kernels[:, i].max(axis=0).toarray() + np.finfo(float).eps)
            k_norm = sparse.hstack([k_norm, kernels[:, i].multiply(inv_max)])
        inv_max = k_norm.max()
        k_norm = k_norm.multiply(inv_max)
        # add a row of ones for the regression
        k_norm = sparse.hstack([k_norm, np.ones((k_norm.shape[0], 1))])

        self.tell('Building combination target...')
        s = self.seed.tocsc()
        s = s/s.max()
        binarize(s, threshold=0.2, copy=False)
        idx = np.asarray(s.sum(axis=0) > 0)
        t = s[:, idx[0]]

        self.tell('Computing Jaccard similarity...')
        j = sparse.coo_matrix(Utilities.jaccard(t))
        jac = Graph.to_sparse_vector(sparse.triu(j, 1)).toarray()

        self.tell('Combining graphs...')
        res = sparse.linalg.lsqr(k_norm, jac)
        self.combined = sparse.coo_matrix(k_norm.tocsc()[:, :-1].dot(res[0][:-1]))
        self.tell('Coefficients:', res[0][:-1])



