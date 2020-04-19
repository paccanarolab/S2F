import os

import pandas as pd
from scipy import sparse

from diffusion.ConsistencyMethod import ConsistencyMethod
from diffusion.LabelWeightedPropagation import LabelWeightedPropagation
from diffusion.S2FLabelPropagation import S2FLabelPropagation
from Utils import ColourClass, FancyApp, Utilities


class Diffuse(FancyApp.FancyApp):

    def __init__(self, args):
        super(Diffuse, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_GREEN
        # this file requires no configuration, only paths
        self.graph_file = args.graph
        self.graph_sep = args.graph_sep
        self.labelling_file = args.labelling
        self.labelling_sep = args.labelling_sep
        self.diffusion_method = args.diffusion_method
        self.kernel_parameters = args.kernel_parameters
        self.graph = None
        self.labelling = None
        self.fasta = args.fasta
        self.output = args.output
        self.proteins = Utilities.extract_indices_from_fasta(
                                self.fasta, Utilities.keep_uniprot_accession)
        self.terms = None

    def run(self):
        self.tell('Selecting', self.diffusion_method, 'as diffusion method')
        self.diffusion_method = self.select_diffusion_method()
        self.tell('Parsing graph...')
        self.read_graph()
        self.tell('Parsing labelling...')
        self.read_labelling()
        diffusion = self.diffusion_method(self.graph, self.proteins,
                                          self.terms)
        diffusion.compute_kernel(**self.kernel_parameters)
        diffusion.diffuse(self.labelling, **self.kernel_parameters)
        self.tell('Saving diffusion results into:', self.output)
        diffusion.write_results(self.output)

    def select_diffusion_method(self):
        if self.diffusion_method == 's2f':
            return S2FLabelPropagation
        elif self.diffusion_method == 'consistency-method':
            return ConsistencyMethod
        elif self.diffusion_method == 'label-weighted':
            return LabelWeightedPropagation
        raise ModuleNotFoundError("Could not find selected " +
                                  "propagation method " +
                                  self.diffusion_method)

    def read_graph(self):
        graph_df = pd.read_csv(os.path.expanduser(self.graph_file),
                               sep=self.graph_sep,
                               names=['protein1', 'protein2', 'score'])
        graph_df = graph_df.merge(self.proteins, left_on='protein1',
                                  right_index=True)
        graph_df = graph_df.merge(self.proteins, left_on='protein2',
                                  right_index=True, suffixes=['1', '2'])
        p1_idx = graph_df['protein idx1'].values
        p2_idx = graph_df['protein idx2'].values

        A = sparse.coo_matrix((graph_df['score'], (p1_idx, p2_idx)),
                              shape=(len(self.proteins),
                                     len(self.proteins))).tocsr()
        self.graph = A.T + A
        # this avoids the efficiency warning, although it can't be too bad
        self.graph = self.graph.tolil()
        self.graph.setdiag(0)  # avoid self-loops
        self.graph = self.graph.tocsc()
        self.tell('Graph dimensions:', self.graph.shape)

    def read_labelling(self):
        labelling_df = pd.read_csv(os.path.expanduser(self.labelling_file),
                                   sep=self.labelling_sep,
                                   names=['protein', 'goterm', 'score'])
        self.terms = pd.DataFrame(list(enumerate(sorted(
                                    labelling_df['goterm'].unique()))))
        self.terms.columns = ['term idx', 'term id']
        self.terms.set_index('term id', inplace=True)
        labelling_df = labelling_df.merge(self.proteins, left_on='protein',
                                          right_index=True)
        labelling_df = labelling_df.merge(self.terms, left_on='goterm',
                                          right_index=True)
        p_idx = labelling_df['protein idx'].values
        go_idx = labelling_df['term idx'].values
        self.labelling = sparse.coo_matrix((labelling_df['score'],
                                            (p_idx, go_idx)),
                                           shape=(len(self.proteins),
                                                  len(self.terms)))
        self.tell('Labelling dimensions:', self.labelling.shape)

