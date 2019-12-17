import os

import pandas as pd

from GeneOntology import GeneOntology

__author__ = 'Mateo Torres'
__email__ = 'Mateo.Torres.2015@live.rhul.ac.uk'
__copyright__ = 'Copyright (c) 2020, Mateo Torres'
__license__ = 'MIT'
__version__ = '0.1'

go = GeneOntology('example data/go.obo')
go.build_structure()
org = 'Acetoanaerobium sticklandii'

organisms = []
data = []
for f in os.listdir('example data'):
    if f.endswith('.goa'):
        d = {}
        organisms.append(os.path.basename(f))
        go.load_annotation_file('example data/'+f, organisms[-1],
                                GeneOntology.ALL_EVIDENCE_CODES)
        go.load_annotation_file('example data/'+f, organisms[-1] + ' exp',
                                GeneOntology.EXPERIMENTAL_EVIDENCE_CODES)
        go.up_propagate_annotations(organisms[-1] + ' exp')

        # any evidence code
        annotated_terms = [t for t in go.terms.values()
                           if organisms[-1] in t.annotations]
        annotated_genes = set()
        for t in annotated_terms:
            annotated_genes |= t.annotations[organisms[-1]]

        d['organism'] = organisms[-1]
        d['all genes'] = len(annotated_genes)

        # experimental annotations
        organism = organisms[-1] + ' exp'
        annotated_terms = [t for t in go.terms.values()
                           if organism in t.annotations]
        annotated_genes = set()
        annotations_by_domain = {'biological_process': set(),
                                 'cellular_component': set(),
                                 'molecular_function': set()}
        terms_by_domain = {'biological_process': set(),
                           'cellular_component': set(),
                           'molecular_function': set()}
        for t in annotated_terms:
            annotated_genes |= t.annotations[organism]
            annotations_by_domain[t.domain] |= t.annotations[organism]
            terms_by_domain[t.domain].add(t)

        # terms annotated with 3 genes or more
        popular_terms = [t for t in annotated_terms
                         if len(t.annotations[organism]) >= 3]
        # genes annotated to those terms
        popular_genes = set()
        popular_by_domain = {'biological_process': set(),
                             'cellular_component': set(),
                             'molecular_function': set()}
        popular_terms_by_domain = {'biological_process': set(),
                                   'cellular_component': set(),
                                   'molecular_function': set()}
        for t in popular_terms:
            popular_genes |= t.annotations[organism]
            popular_by_domain[t.domain] |= t.annotations[organism]
            popular_terms_by_domain[t.domain].add(t)

        d['exp annotations'] = len(annotated_genes)
        d['popular genes'] = len(popular_genes)
        d['annotated terms'] = len(annotated_terms)
        d['popular terms'] = len(popular_terms)
        d['terms bp'] = len(terms_by_domain['biological_process'])
        d['terms mf'] = len(terms_by_domain['molecular_function'])
        d['terms cc'] = len(terms_by_domain['cellular_component'])
        d['pop bp'] = len(popular_terms_by_domain['biological_process'])
        d['pop mf'] = len(popular_terms_by_domain['molecular_function'])
        d['pop cc'] = len(popular_terms_by_domain['cellular_component'])

        data.append(d)

df = pd.DataFrame(data)
df.to_pickle('data.pck')
# 1956,6,6,66,20,40,26,0,18,2,0
