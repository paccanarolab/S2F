from Utils import *
from graphs import Graph
from scipy import sparse

import pandas as pd
import numpy as np
import os
import subprocess


class Homology(Graph):

    def __init__(self, fasta, proteins, graphs_dir, alias):
        super(Homology, self).__init__()
        self.fasta = fasta
        self.proteins = proteins
        self.graphs_dir = graphs_dir
        self.alias = alias
        self.homology_graph = None

    def get_graph(self, **kwargs):
        h = self.homology_graph.merge(self.proteins, left_on='Protein 1', right_index=True)
        h = h.merge(self.proteins, left_on='Protein 2', right_index=True, suffixes=['1', '2'])
        p1_idx = h['protein idx1'].values
        p2_idx = h['protein idx2'].values
        return sparse.coo_matrix((h['weight'], (p1_idx, p2_idx)),
                                 shape=(len(self.proteins), len(self.proteins)))

    def write_graph(self, filename):
        Graph.assert_lexicographical_order(self.homology_graph)
        self.homology_graph.to_csv(filename, sep='\t')

    def compute_graph(self):
        homology_graph = os.path.join(self.graphs_dir, self.alias)
        if not os.path.exists(homology_graph):
            self.tell('Computing homology graph...')
            # compute the homology graph
            blast_command = 'blastp -query {fasta} -db {blastdb} -out {out} -outfmt "6 std qlen"'
            out = os.path.join(self.graphs_dir, self.alias + '_homology.blast')
            if not os.path.exists(out):
                self.tell(blast_command.format(fasta=self.fasta, blastdb=self.fasta, out=out))
                subprocess.call(blast_command.format(fasta=self.fasta, blastdb=self.fasta, out=out), shell=True)

            self.tell('Parsing BLAST output...')
            # parse the blast output
            homology = {}
            for line in open(out):
                fields = line.strip().split('\t')
                query_id = fields[0]
                subject_id = fields[1]
                evalue = float(fields[10])

                if query_id not in homology:
                    homology[query_id] = {}

                if subject_id in homology[query_id]:
                    homology[query_id][subject_id] = np.min([homology[query_id][subject_id], evalue])
                else:
                    homology[query_id][subject_id] = evalue

            self.tell('Building graph...')
            graph = {}
            hom_graph = {'Protein 1': [], 'Protein 2': [], 'weight': []}
            ind_matrix = []
            proteins = sorted(homology.keys())
            maxi = -1.0
            for i, p1 in enumerate(proteins):
                for p2 in proteins[i:]:

                    e12 = 10.0
                    e21 = 10.0

                    if p2 in homology[p1]:
                        e12 = homology[p1][p2]
                    if p1 in homology[p2]:
                        e21 = homology[p2][p1]

                    val = np.max([e12, e21])
                    transformed = 0.0

                    if val == 0 or p1 == p2:
                        ind_matrix.append((p1, p2))
                        ind_matrix.append((p2, p1))
                        transformed = 1.0
                    elif 0.0 < val < 11.0:
                        transformed = -1.0 * np.log(val/11.0)
                        maxi = np.max([maxi, transformed])
                    graph[p1, p2] = transformed
                    graph[p2, p1] = transformed
            inv_max = 1.0 / maxi

            self.tell('Normalising homology graph...')
            for i, p1 in enumerate(proteins):
                for p2 in proteins[i:]:
                    val = graph[p1, p2]
                    hom_graph['Protein 1'].append(p1)
                    hom_graph['Protein 2'].append(p2)
                    if (p1, p2) in ind_matrix:
                        hom_graph['weight'].append(val)
                    else:
                        hom_graph['weight'].append(inv_max * val)
            self.homology_graph = pd.DataFrame.from_dict(hom_graph)

            # we make sure that lexicographical order is maintained so that
            # all values are kept on the upper triangle in the matrix
            Graph.assert_lexicographical_order(self.homology_graph)

            self.homology_graph.to_pickle(homology_graph)
        else:
            self.tell('Homology graph found, skipping computation...')
            self.homology_graph = pd.read_pickle(homology_graph)
            Graph.assert_lexicographical_order(self.homology_graph)





