import os

import numpy as np
import pandas as pd
from scipy import sparse

from seeds import Seed
from Utils import ColourClass, Utilities
from Utils.HmmerTbloutParser import HmmerTbloutFile


class HMMerSeed(Seed):

    def __init__(self, hmmer, proteins, terms, go, blacklist, goa, protein_format):
        super(HMMerSeed, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_CYAN
        self.hmmer = hmmer
        self.proteins = proteins
        self.terms = terms
        self.go = go
        self.blacklist = blacklist
        self.goa = goa
        self.evalue_file = ''
        self.seed_threshhold = 1e-4
        self.protein_format = protein_format

    def get_seed(self, **kwargs):
        if 'seed_threshhold' in kwargs:
            self.seed_threshhold = kwargs['seed_threshhold']
        return_pandas_assignment = kwargs.get('return_pandas_assignment',
                                              False)
        assignment = {'GO ID': [], 'Protein': [], 'Score': []}
        self.tell('Building HMMer seed file with threshold',
                  self.seed_threshhold)
        for line in open(self.evalue_file, 'r'):
            fields = line.strip().split('\t')
            evalue = float(fields[2])
            if evalue <= self.seed_threshhold and len(fields) >= 4:
                prot = fields[0]
                for t in fields[3:]:
                    assignment['GO ID'].append(t)
                    assignment['Protein'].append(prot)
                    assignment['Score'].append(1)
        assignment_df = pd.DataFrame(assignment)
        self.tell('Up-propagating assignment...')
        self.go.load_annotations(assignment_df, 'HMMer seed')
        self.go.up_propagate_annotations('HMMer seed')
        assignment_df = self.go.get_annotations('HMMer seed')

        if return_pandas_assignment:
            return assignment_df

        self.tell('Building seed matrix...')
        assignment_df = assignment_df.merge(self.proteins, left_on='Protein',
                                            right_index=True)
        assignment_df = assignment_df.merge(self.terms, left_on='GO ID',
                                            right_index=True)

        p_idx = assignment_df['protein idx'].values
        go_idx = assignment_df['term idx'].values

        return sparse.coo_matrix((np.ones(len(assignment_df)),
                                  (p_idx, go_idx)),
                                 shape=(len(self.proteins),
                                        len(self.terms)))

    def process_output(self, **kwargs):
        self.evalue_file = kwargs.get('evalue_file', self.evalue_file)
        if not os.path.exists(self.evalue_file):
            self.tell('Loading GOA annotations')
            self.go.load_annotation_file(self.goa, 'GOA',
                                         blacklist=self.blacklist)
            self.tell('Caching GOA annotations')
            goa = self.go.get_annotations('GOA')
            self.tell('Parsing HMMER file')
            hmmer_parser = HmmerTbloutFile(self.hmmer)
            data = {k:[] for k in ['target', 'query', 'evalue']}
            for h in hmmer_parser:
                target = h.target_name
                query = h.query_name
                if self.protein_format == 'uniprot':
                    query = Utilities.extract_uniprot_accession(query) 
                target = Utilities.extract_uniprot_accession(target)
                data['target'].append(target)
                data['query'].append(query)
                data['evalue'].append(h.e_value)
            hmmer_df = pd.DataFrame(data)
            del data
            self.tell('Transferring function from GOA')
            m = hmmer_df.merge(goa[['Protein', 'GO ID']], 
                               left_on='target', 
                               right_on='Protein')[['target', 
                                                    'query', 
                                                    'evalue', 
                                                    'GO ID']]
            m = m.groupby(['query', 'target', 'evalue'], as_index = False).agg({'GO ID': '\t'.join})
            # the line below can be a bit slower than aggregating a list of GO terms in pandas
            # but it generates way fewer empty columns
            self.tell('Preparing evalue file')
            m['write'] = m.apply(lambda r: '\t'.join(r.dropna().astype(str).values), axis=1)
            self.tell('Writing evalue file')
            m['write'].to_csv(self.evalue_file, index=False, header=False)
            
        else:
            self.tell('evalue file found, skipping computation')
