#!/usr/bin/env python
"""

An utility Script and tool for managing Gene Ontology related
work.

=======
License
=======

Copyright (c) 2018 Mateo Torres <Mateo.Torres.2016@live.rhul.ac.uk>

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
"""


from collections import defaultdict

import numpy as np
import pandas as pd

from Utils import ColourClass, FancyApp, ProgressBar

if __name__ == '__main__':
    import OboParser
    import AnnotationParser
else:
    from . import OboParser
    from . import AnnotationParser

__author__ = 'Mateo Torres'
__email__ = 'Mateo.Torres.2015@live.rhul.ac.uk'
__copyright__ = 'Copyright (c) 2018, Mateo Torres'
__license__ = 'MIT'
__version__ = '0.1'


class GOTerm(object):

    VALID_RELATIONS = ['is_a', 'part_of']

    def __init__(self, go_id, name, domain, ontology):
        self.go_id = go_id
        self.name = name
        self.relations = defaultdict(set)
        self.annotations = defaultdict(dict)
        self.domain = domain
        self.aliases = set()
        self.is_obsolete = False
        self.ic = {}
        self.ontology = ontology

    def information_content(self, organism_name):
        if organism_name not in self.ic:
            annotations_df = self.ontology.get_annotations(organism_name)
            len_annotations = len(self.annotations[organism_name])
            if len_annotations > 0:
                self.ic[organism_name] = -np.log(len_annotations /
                                                 annotations_df.shape[0]) /\
                                         np.log(2)
            else:
                self.ic[organism_name] = 0
        return self.ic[organism_name]

    def add_relation(self, go_term, relation_type):
        if relation_type in GOTerm.VALID_RELATIONS:
            self.relations[relation_type].add(go_term)
            go_term.relations['a_'+relation_type].add(self)

    def up_propagate_annotations(self, organism_name,
                                 relations=VALID_RELATIONS,
                                 same_domain=True):
        """
        Recursively up-propagates the annotations until the root term

        Parameters
        ----------
        organism_name : str
            The annotation set to up-propagate
        relations : list, optional
            a list of relations that will be considered during up-propagation,
            (defaults to ['is_a', 'part_of']).
            All relations are assumed to be transitive
        same_domain : bool, optional
            If True, the up-propagation is constrained to terms belonging
            to the same domain
        """
        for relation in relations:
            for term in self.relations[relation]:
                if (same_domain and term.domain == self.domain) or\
                        not same_domain:
                    for prot, score in self.annotations[organism_name].items():
                        if prot in term.annotations[organism_name].keys():
                            term.annotations[organism_name][prot] = max(
                                term.annotations[organism_name][prot],
                                score)
                        else:
                            term.annotations[organism_name][prot] = score
                    term.up_propagate_annotations(organism_name,
                                                  relations=relations,
                                                  same_domain=same_domain)

    def get_parents(self, relations=VALID_RELATIONS):
        parents = set()
        for relation in relations:
            parents |= self.relations[relation]
        return parents

    def get_ancestors(self, relations=VALID_RELATIONS):
        ancestors = set()
        for relation in relations:
            ancestors |= self.relations[relation]
            for term in self.relations[relation]:
                ancestors |= term.get_ancestors(relations)
        return ancestors


class GeneOntology(FancyApp.FancyApp):

    EXPERIMENTAL_EVIDENCE_CODES = ['EXP', 'IDA', 'IPI', 'IMP',
                                   'IGI', 'IEP', 'TAS', 'IC']
    ALL_EVIDENCE_CODES = [
        # experimental
        'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP',
        # High Throughput
        'HTP', 'HDA', 'HMP', 'HGI', 'HEP',
        # Computational Analysis
        'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA'
        # Author statement
        'TAS', 'NAS',
        # Curator statement
        'IC', 'ND',
        # Electronic Annotation
        'IEA'
    ]
    DOMAINS = ['biological_process',
               'cellular_component',
               'molecular_function']
    RELATIONS = ['is_a', 'part_of']

    def __init__(self, obo, verbose=False):
        super(GeneOntology, self).__init__()
        self.colour = ColourClass.bcolors.OKGREEN
        self.obo = open(obo)
        self.terms = {}
        self.alias_map = {}  # a cache of the aliases for faster access
        self.__verbose__ = verbose

        # Some useful references

    def find_term(self, go_id):
        """
        If the go_id is in the structure, return the term, otherwise,
        find by alias
        :param go_id:
        :return: GOTerm instance
        """
        try:
            return self.terms[go_id]
        except KeyError:
            return self.terms[self.alias_map[go_id]]

    def build_structure(self):
        """
        Builds the structure of the Gene Ontology in memory
        based on the provided obo file.
        :return: None
        """
        parser = OboParser.Parser(self.obo)

        stanzas = []
        self.tell('parsing Obo File')
        for stanza in parser:
            # Ignore every stanza not describing a Term
            if stanza.name != 'Term':
                continue

            # Create the term
            term = GOTerm(
                stanza.tags['id'][0].value,
                stanza.tags['name'][0].value,
                stanza.tags['namespace'][0].value,
                self
            )
            self.terms[term.go_id] = term

            # Set the aliases..
            if 'alt_id' in stanza.tags:
                for alias in stanza.tags['alt_id']:
                    self.alias_map[alias.value] = term.go_id
                    self.terms[term.go_id].aliases.add(alias.value)

            # set the 'is obsolete' variable
            if 'is_obsolete' in stanza.tags:
                self.terms[term.go_id].is_obsolete = True

            stanzas.append(stanza)

        self.tell('Building Structure')
        if self.__verbose__:
            prog = ProgressBar.ProgressBar(0, len(stanzas), 77, mode='dynamic',
                                           char='-')
        for stanza in stanzas:
            go_id = stanza.tags['id'][0].value
            if 'is_a' in stanza.tags:
                for related_go_id in stanza.tags['is_a']:
                    self.find_term(go_id)\
                        .add_relation(self.find_term(related_go_id.value),
                                      'is_a')

            if 'relationship' in stanza.tags:
                for relationship in stanza.tags['relationship']:
                    if relationship.modifiers is not None:
                        self.warning('modifiers is not None')
                    split_relation = relationship.value.split()
                    if split_relation[0] == 'part_of':
                        self.find_term(go_id)\
                            .add_relation(self.find_term(split_relation[1]),
                                          'part_of')
            if self.__verbose__:
                prog.increment_amount()
                prog.print_bar()
        if self.__verbose__:
            prog.finish_bar()

    def load_clamp_file(self, clamp_file, organism_name, annotate_obsolete=False):
        """
        Loads a tab separated annotation file into the structure. The structure
        of such file mist be:
        <protein_id>TAB<go_term_id>
        Parameters
        ----------
        clamp_file : path to the file
        organism_name : the name that keeps these annotations
                              separated from other annotations in the structure
        annotate_obsolete : bool, optional
            if True, obsolete terms will be considered in the annotation
            (the default is False)
        """
        for line in open(clamp_file):
            fields = line.strip().split('\t')
            term = self.find_term(fields[1])
            if annotate_obsolete or not term.is_obsolete:
                term.annotations[organism_name][fields[0]] = 1

    def load_annotation_file(self, annotation_file, organism_name,
                             evidence_codes=EXPERIMENTAL_EVIDENCE_CODES,
                             domains=DOMAINS, blacklist=None,
                             annotate_obsolete=False):
        """
        Loads a GOA annotation file into the structure.

        Parameters
        ----------
        annotation_file : path to the GOA file
        organism_name : the name that keeps these annotations
                              separated from other annotations in the structure
        evidence_codes : a list of valid evidence codes, only
                               anntoations with evidence codes in the list
                               will be added to the structure. Defaults to
                               `EXPERIMENTAL_EVIDENCE_CODES`
        domains : a list of domains that will be considered, defaults to
                        `DOMAINS`
        blacklist : list of taxons to ignore. Annotations associated to
                          these taxons will be gnored.
        annotate_obsolete : bool, optional
            if True, obsolete terms will be considered in the annotation
            (the default is False)
        """
        parser = AnnotationParser.AnnotationFile(annotation_file,
                                                 organism_name)
        counter = 0
        self.tell('Annotating ' + organism_name)
        for annotation in parser:
            # check that the organism is not blacklisted
            if blacklist:
                blacklisted = [i for i in annotation.taxons if i in blacklist]
                if blacklisted:
                    continue
            # check the desired evidence code
            if annotation.evidence_code in evidence_codes:
                try:
                    term = self.find_term(annotation.go_id)
                    if term.domain in domains:
                        if len(annotation.qualifiers) > 0:
                            if annotation.qualifiers[0] != 'NOT':
                                if annotate_obsolete or not term.is_obsolete:
                                    term.annotations[organism_name][
                                        annotation.db_object_id] = 1
                            else:
                                counter += 1
                        else:
                            if annotate_obsolete or not term.is_obsolete:
                                term.annotations[organism_name][
                                    annotation.db_object_id] = 1
                    else:
                        counter += 1
                except KeyError:
                    counter += 1

        self.warning('A total of ' + str(counter) +
                     ' annotations were skipped')

    def load_annotations(self, annotations, organism_name,
                         annotate_obsoletes=False):
        """
        Load annotations from a compatible pandas DataFrame

        Parameters
        ----------
        annotations : pandas DataFrame
            with the columns 'GO ID', 'Protein' and 'Score'
        organism_name : str
            the name of the organism in which to load the annotations,
            if the organism exists, annotations will be overwritten
        annotate_obsoletes : bool, optional
            if True, obsolete terms will be considered in the annotation
            (the default is False)

        Notes
        -----
        The expected DataFrame should contain only one score per pair,
        otherwise the latest score will be kept.
        """
        for i, a in annotations.iterrows():
            term = self.find_term(a['GO ID'])
            if annotate_obsoletes or not term.is_obsolete:
                term.annotations[organism_name][a['Protein']] = a['Score']

    def up_propagate_annotations(self, organism_name, relations=RELATIONS):
        annotated_terms = set()
        # retrieve annotated terms
        for term in self.terms.values():
            if organism_name in term.annotations:
                annotated_terms.add(term)

        # up propagate
        for term in annotated_terms:
            term.up_propagate_annotations(organism_name, relations=relations)

    def dump_annotations(self, organism_name, out, format='pandas'):
        if format == 'pandas':
            self.get_annotations(organism_name).to_pickle(out)
        elif format == 'txt':
            outfile = open(out)
            for term in self.terms.values():
                if organism_name in term.annotations:
                    for protein in term.annotations[organism_name]:
                        outfile.write(term.go_id + '\t' + protein + '\n')
            outfile.close()

    def get_annotations(self, organism_name):
        d = {'GO ID': [], 'Protein': [], 'Score': []}
        for term in self.terms.values():
            if organism_name in term.annotations:
                for protein, score in term.annotations[organism_name].items():
                    d['GO ID'].append(term.go_id)
                    d['Protein'].append(protein)
                    d['Score'].append(score)
        return pd.DataFrame(d)


if __name__ == '__main__':
    go = GeneOntology('example data/go.obo', verbose=True)
    go.build_structure()
    go.load_annotation_file('example data/36.P_aeruginosa_LMG_12228.goa',
                            'Pseudomonas aeruginosa')

    annotated_terms = [t for t in go.terms.values() if 'Pseudomonas aeruginosa'
                       in t.annotations.keys()]
    go.tell('annotated terms: ', len(annotated_terms))
    go.tell('number of annotations: ',
            sum([len(t.annotations['Pseudomonas aeruginosa'])
                 for t in annotated_terms]))

    go.dump_annotations('Pseudomonas aeruginosa',
                        'example data/pseudomonas.exp')
    go.up_propagate_annotations('Pseudomonas aeruginosa')
    go.dump_annotations('Pseudomonas aeruginosa',
                        'example data/pseudomonas_uppropagated.exp')

    annotated_terms = [t for t in go.terms.values()
                       if 'Pseudomonas aeruginosa'
                       in t.annotations.keys()]
    go.tell('number of annotations: ',
            sum([len(t.annotations['Pseudomonas aeruginosa'])
                 for t in annotated_terms]))

    annotated_genes = set()
    terms_by_domain = {'biological_process': set(),
                       'cellular_component': set(),
                       'molecular_function': set()}
    for t in annotated_terms:
        annotated_genes |= set(t.annotations['Pseudomonas aeruginosa'].keys())
        terms_by_domain[t.domain].add(t)
    go.tell('number of annotated genes: ', len(annotated_genes))

    # terms annotated with 3 genes or more
    popular_terms = [t for t in annotated_terms
                     if len(t.annotations['Pseudomonas aeruginosa'].keys())
                     >= 3]
    # genes annotated to those terms
    popular_genes = set()
    popular_by_domain = {'biological_process': set(),
                         'cellular_component': set(),
                         'molecular_function': set()}
    popular_terms_by_domain = {'biological_process': set(),
                               'cellular_component': set(),
                               'molecular_function': set()}
    for t in popular_terms:
        popular_genes |= set(t.annotations['Pseudomonas aeruginosa'].keys())
        popular_by_domain[t.domain] |= set(
            t.annotations['Pseudomonas aeruginosa'].keys())
        popular_terms_by_domain[t.domain].add(t)
    go.tell('popular genes', len(popular_genes))
    go.tell('annotated terms: ', len(annotated_terms))
    go.tell('popular terms: ', len(popular_terms))

    go.tell('terms bp', len(terms_by_domain['biological_process']))
    go.tell('terms mf', len(terms_by_domain['molecular_function']))
    go.tell('terms cc', len(terms_by_domain['cellular_component']))

    go.tell('pop terms bp', len(popular_terms_by_domain['biological_process']))
    go.tell('pop terms mf', len(popular_terms_by_domain['molecular_function']))
    go.tell('pop terms cc', len(popular_terms_by_domain['cellular_component']))

    go.tell('popular genes bp', len(popular_by_domain['biological_process']))
    go.tell('popular genes mf', len(popular_by_domain['molecular_function']))
    go.tell('popular genes cc', len(popular_by_domain['cellular_component']))
