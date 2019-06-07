from Utils import *
from graphs import Graph
from graphs import combination

from scipy import sparse
from scipy.sparse.linalg import lsqr
from sklearn.preprocessing import binarize

import pandas as pd
import numpy as np
import os


class Combine(FancyApp.FancyApp):

    def __init__(self, args):
        super(Combine, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_GREEN
        # this file requires no configuration, only paths
        self.collection_file = args.collection
        self.collection_sep = args.collection_sep
        self.collection_selection = args.collection_selection
        self.homology_file = args.homology
        self.homology_sep = args.homology_sep
        self.seed_file = args.seed
        self.seed_sep = args.seed_sep
        self.seed_threshold = args.seed_threshold
        self.fasta = args.fasta
        self.output = args.output
        self.collection = {}
        self.homology = None
        self.seed = None
        self.combined = None
        self.proteins = Utilities.extract_indices_from_fasta(self.fasta, Utilities.keep_uniprot_accession)

    def run(self):
        self.tell('parsing graph combination file...')
        self.collection = self.read_graph_collection()
        self.tell('parsing homology graph file...')
        self.homology = self.read_homology_file()
        self.tell('parsing seed file...')
        self.seed = self.get_seed_similarity()
        comb = combination.Combination(self.proteins, self.collection, self.homology, self.seed)
        comb.compute_graph()
        self.tell('writing output file:', self.output)
        comb.write_graph(self.output)
        self.tell('done')

    def read_graph_collection(self):
        # 1. read collection file
        collection_df = pd.read_csv(self.collection_file, sep=self.collection_sep)
        cols = collection_df.columns

        # 2. manually set score columns to float
        collection_df[cols[2:]] = collection_df[cols[2:]].astype('float32')

        # 3. we guarantee the lexicographical order between the protein columns
        Graph.assert_lexicographical_order(collection_df, p1=cols[0], p2=cols[1])

        # 4. forget any self-similarity
        collection_df = collection_df[collection_df[cols[0]] != collection_df[cols[1]]]

        # 5. build the graph collection
        collection = {}
        for g in self.collection_selection:
            self.tell('processing', g)
            c = list(cols[:2])
            c.append(g)
            graph = collection_df[c]
            graph = graph[graph[g]>0]
            Graph.assert_lexicographical_order(graph, cols[0], cols[1])

            graph = graph.merge(self.proteins, left_on=cols[0], right_index=True)
            graph = graph.merge(self.proteins, left_on=cols[1], right_index=True, suffixes = ['1', '2'])

            p1_idx = graph['protein idx1'].values
            p2_idx = graph['protein idx2'].values

            collection[g] = sparse.coo_matrix((graph[g], (p1_idx, p2_idx)),
                                              shape=(len(self.proteins), len(self.proteins)))
        return collection

    def get_seed_similarity(self):
        df = pd.read_csv(self.seed_file, sep=self.seed_sep,
                         names=["Protein", "GO ID", "Score"],
                         dtype={"Score": 'float32',
                                "Protein": "str",
                                "GO ID": "str"})
        df = df[df['Score'] > 0]
        df = df.merge(self.proteins, left_on='Protein', right_index=True)
        t = pd.DataFrame(list(enumerate(np.sort(df['GO ID'].unique()))))
        t.columns = ['term idx', 'term id']
        t.set_index('term id', inplace=True)
        df = df.merge(t, left_on='GO ID', right_index=True)
        p_idx = df['protein idx'].values
        go_idx = df['term idx'].values
        return sparse.coo_matrix((df['Score'].values, (p_idx, go_idx)),
                                 shape=(len(self.proteins), len(t)))

    def read_homology_file(self):
        columns = ['protein1', 'protein2', 'score']
        types = dict([(x, 'float32') if "protein" not in x
                      else (x, "str") for x in columns])
        df = pd.read_csv(self.homology_file, sep=self.homology_sep, names=columns, dtype=types)

        min_protein = df[["protein1", "protein2"]].min(axis=1)
        max_protein = df[["protein1", "protein2"]].max(axis=1)
        df.loc[:, "protein1"] = min_protein
        df.loc[:, "protein2"] = max_protein

        df = df.merge(self.proteins, left_on='protein1', right_index=True)
        df = df.merge(self.proteins, left_on='protein2', right_index=True, suffixes=['1', '2'])
        p1_idx = df['protein idx1'].values
        p2_idx = df['protein idx2'].values
        return sparse.coo_matrix((df['score'], (p1_idx, p2_idx)),
                                 shape=(len(self.proteins), len(self.proteins)))