from Utils import *
from seeds import Seed
from scipy import sparse

import pandas as pd
import numpy as np


class InterProSeed(Seed):

    def __init__(self, interpro, proteins, terms, go):
        super(InterProSeed, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_CYAN
        self.interpro = interpro
        self.proteins = proteins
        self.terms = terms
        self.go = go
        self.all_methods = None
        self.methods = []

    def get_seed(self):
        self.tell('Combining', len(self.methods), 'models in InterPro')
        return self.all_methods * (1.0/len(self.methods))

    def process_output(self, **kwargs):
        methods = {}
        self.tell('Processing InterPro output')
        for line in open(self.interpro):
            fields = line.split('\t')
            method = fields[3].lower()
            prot = fields[0]
            if len(fields) >= 14:
                if method not in ['seg', 'coil']:
                    terms = fields[13].strip()
                    if len(terms) > 0:
                        terms = terms.split('|')
                        for t in terms:
                            if method in methods.keys():
                                methods[method]['GO ID'].append(t)
                                methods[method]['Protein'].append(prot)
                                methods[method]['Score'].append(1)
                            else:
                                self.methods.append(method)
                                methods[method] = {'GO ID': [t], 'Protein': [prot], 'Score': [1]}

        self.methods.sort()
        methods_df = {}
        self.all_methods = sparse.csr_matrix((len(self.proteins), len(self.terms)))
        self.tell('Annotating and up propagating')
        for e, m in enumerate(self.methods):
            methods_df[m] = pd.DataFrame(methods[m])
            self.go.load_annotations(methods_df[m], m)
            self.go.up_propagate_annotations(m)
            methods_df[m] = self.go.get_annotations(m)

            methods_df[m] = methods_df[m].merge(self.proteins, left_on='Protein', right_index=True)
            methods_df[m] = methods_df[m].merge(self.terms, left_on='GO ID', right_index=True)

            p_idx = methods_df[m]['protein idx'].values
            go_idx = methods_df[m]['term idx'].values

            self.all_methods += sparse.coo_matrix((np.ones(len(methods_df[m])), (p_idx, go_idx)),
                                                  shape=(len(self.proteins), len(self.terms)))

        return methods_df


