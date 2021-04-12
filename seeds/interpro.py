import numpy as np
import pandas as pd
from scipy import sparse
from seeds import Seed
from Utils import ColourClass, Utilities


class InterProSeed(Seed):

    def __init__(self, interpro, proteins, terms, go, protein_format):
        super(InterProSeed, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_CYAN
        self.interpro = interpro
        self.proteins = proteins
        self.terms = terms
        self.go = go
        self.protein_format = protein_format
        self.all_methods = None
        self.methods = []

    def get_seed(self, **kwargs):
        return_pandas_assignment = kwargs.get('return_pandas_assignment',
                                              False)
        self.tell('Combining', len(self.methods), 'models in InterPro')
        combined = self.all_methods * (1.0/len(self.methods))

        if return_pandas_assignment:
            return Seed._seed_to_pandas(combined.tocoo(),
                                        self.proteins,
                                        self.terms)
        return combined

    def process_output(self, **kwargs):
        methods = {}
        self.tell('Processing InterPro output')
        for line in open(self.interpro):
            fields = line.split('\t')
            method = fields[3].lower()
            prot = fields[0]
            if self.protein_format is 'uniprot':
                prot = Utilities.extract_uniprot_accession(prot) 
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
                                methods[method] = {'GO ID': [t],
                                                   'Protein': [prot],
                                                   'Score': [1]}

        self.methods.sort()
        methods_df = {}
        self.all_methods = sparse.csr_matrix((len(self.proteins),
                                              len(self.terms)))
        self.tell('Annotating and up propagating')
        for e, m in enumerate(self.methods):
            methods_df[m] = pd.DataFrame(methods[m])
            self.go.load_annotations(methods_df[m], m)
            self.go.up_propagate_annotations(m)
            methods_df[m] = self.go.get_annotations(m)

            methods_df[m] = methods_df[m].merge(self.proteins,
                                                left_on='Protein',
                                                right_index=True)
            methods_df[m] = methods_df[m].merge(self.terms, left_on='GO ID',
                                                right_index=True)

            p_idx = methods_df[m]['protein idx'].values
            go_idx = methods_df[m]['term idx'].values

            self.all_methods += sparse.coo_matrix((np.ones(len(methods_df[m])),
                                                   (p_idx, go_idx)),
                                                  shape=(len(self.proteins),
                                                         len(self.terms)))

        return methods_df
