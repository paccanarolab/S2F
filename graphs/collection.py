import os
import subprocess
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import sparse

from graphs import Graph
from Utils import ColourClass, FancyApp, Utilities


def parse_blast(filename, protein_format):
    best_hit_dict = {}
    for line in open(filename):
        parts = line.strip().split('\t')
        query = parts[0]
        if protein_format == 'uniprot':
            query = Utilities.extract_uniprot_accession(query)
        evalue = parts[10]
        if query in best_hit_dict and best_hit_dict[query]['evalue'] < evalue:
            continue

        target = parts[1].strip()
        if protein_format == 'uniprot':
            target = Utilities.extract_uniprot_accession(target)
        pi = float(parts[2].strip())
        perc = 100.0 * float(parts[7]) - float(parts[6]) / float(parts[-1])
        best_hit_dict[query] = {
            'evalue': evalue,
            'pi': pi,
            'perc': perc,
            'target': target
        }
    return best_hit_dict


def compute_ortholog(fasta, db, string_fasta, string_db, out, out_dir, protein_format):
    outfile = os.path.join(out_dir, out)
    if not os.path.exists(outfile):
        blast_command = 'blastp -query {fasta} -db {blastdb} -out {out}' +\
                        ' -outfmt "6 std qlen" -max_target_seqs 1'

        # create temporary blast files
        out_1 = os.path.join(out_dir, out + 'forward')
        out_2 = os.path.join(out_dir, out + 'backward')
        subprocess.call(blast_command.format(fasta=fasta, blastdb=string_db,
                                             out=out_1), shell=True)
        subprocess.call(blast_command.format(fasta=string_fasta, blastdb=db,
                                             out=out_2), shell=True)

        # max_evalue = 1e-6
        # perc_th = 80.0
        # positives = 60.0

        orthologs_found = {}

        forward = parse_blast(out_1, protein_format)
        backward = parse_blast(out_2, protein_format)
        index = 0
        for query, infoq in forward.items():
            if infoq['target'] in backward.keys():
                infot = backward[infoq['target']]
                if infot['target'] == query:  # reciprocal best hit!
                    orthologs_found[index] = {
                        'query': query,
                        'target': infoq['target'],
                        'query_evalue': infoq['evalue'],
                        'query_pi': infoq['pi'],
                        'query_perc': infoq['perc'],
                        'target_evalue': infot['evalue'],
                        'target_pi': infot['pi'],
                        'target_perc': infot['perc'],
                    }
                    index += 1
        orth = pd.DataFrame.from_dict(orthologs_found, orient='index')
        orth['max evalue'] = np.max(orth[['query_evalue', 'target_evalue']],
                                    axis=1)
        orth['pos'] = np.sqrt(orth['query_pi'] * orth['target_pi'])
        orth.to_pickle(os.path.join(out_dir, out))
        # remove temporary files
        os.remove(out_1)
        os.remove(out_2)
    FancyApp.FancyApp.yell(ColourClass.bcolors.BOLD_CYAN, 'compute ortholog',
                           'Finished collecting orthologs:', out)


class Collection(Graph):

    def __init__(self,
                 fasta, proteins, string_dir, string_links, core_ids,
                 output_dir, orthologs_dir, graphs_dir, alias, cpus, blacklist,
                 max_evalue, perc, positives, protein_format,
                 interesting_graphs=['neighborhood', 'experiments',
                                     'coexpression', 'textmining',
                                     'database']):
        super(Collection, self).__init__()
        self.fasta = fasta
        self.proteins = proteins
        self.string_dir = string_dir
        self.string_links = string_links
        self.core_ids = core_ids
        self.output_dir = output_dir
        self.graphs_dir = graphs_dir
        self.orthologs_dir = orthologs_dir
        self.alias = alias
        self.cpus = cpus
        self.db = self.fasta
        self.blacklist = blacklist
        if self.blacklist is None:
            self.blacklist = []
        self.max_evalue = max_evalue
        self.perc = perc
        self.positives = positives
        self.interesting_graphs = interesting_graphs
        self.collection = None
        self.protein_format = protein_format

    def get_graph(self, **kwargs):
        # It is assumed that the collection was already made
        # so we can focus on the creation of the protein-protein graph only.
        collection = {}
        for g in self.interesting_graphs:
            graph = self.collection[['query1', 'query2', g]]\
                        .merge(self.proteins, left_on='query1',
                               right_index=True)
            graph = graph.merge(self.proteins, left_on='query2',
                                right_index=True, suffixes=['1', '2'])
            p1_idx = graph['protein idx1'].values
            p2_idx = graph['protein idx2'].values
            collection[g] = sparse.coo_matrix((graph[g], (p1_idx, p2_idx)),
                                              shape=(len(self.proteins),
                                                     len(self.proteins)))
        return collection

    def should_be_processed(self, string_id):
        # self.tell('checking', string_id)
        if string_id in self.blacklist:
            # self.tell('blacklisted')
            return False
        orthologs_filename = os.path.join(self.orthologs_dir,
                                          self.alias + '_AND_' + string_id)
        if not os.path.exists(orthologs_filename):
            # self.tell('ortholog file does not exist')
            return False
        orthologs = pd.read_pickle(orthologs_filename)
        if orthologs.shape[0] <= 2:
            # self.tell('not enough orthologs')
            return False
        # self.tell(string_id, 'Should be processed')
        return True

    def clean_graph(self):
        graph = {'protein 1': [], 'protein 2': []}
        for g in self.interesting_graphs:
            graph[g] = []
        return graph

    def process_graph(self, string_id, graph):
        self.tell('Transferring links from', string_id)
        # load ortholog file
        orthologs = pd.read_pickle(os.path.join(self.orthologs_dir,
                                                self.alias + '_AND_' +
                                                string_id))
        orthologs['target_evalue'] = orthologs['target_evalue'].astype(float)
        orthologs['query_evalue'] = orthologs['query_evalue'].astype(float)
        orthologs['max_evalue'] = orthologs['query_evalue'].astype(float)

        valid_orthologs = orthologs[
            (orthologs['max_evalue'] < self.max_evalue) &
            (orthologs['target_perc'] >= self.perc) &
            (orthologs['query_perc'] >= self.perc) &
            (orthologs['pos'] >= self.positives)
        ].copy()

        # valid_orthologs['protein 1'] = valid_orthologs['target']
        # valid_orthologs['protein 2'] = valid_orthologs['target']

        valid_edges = graph[
            (graph['protein 1'].isin(valid_orthologs['target'])) &
            (graph['protein 2'].isin(valid_orthologs['target']))
        ]

        valid_edges = valid_edges.merge(valid_orthologs[['query', 'target',
                                                         'max_evalue']],
                                        left_on='protein 1', right_on='target')
        valid_edges = valid_edges.merge(valid_orthologs[['query', 'target',
                                                         'max_evalue']],
                                        left_on='protein 2', right_on='target',
                                        suffixes=['1', '2'])
        valid_edges['max_evalue'] = valid_edges[['max_evalue1',
                                                 'max_evalue2']].max(axis=1)

        self.tell(valid_edges.shape[0], 'edges transferred from organism',
                  string_id)

        drop_cols = ['protein 1', 'protein 2', 'target1', 'target2',
                     'max_evalue1', 'max_evalue2']
        if self.collection is None:
            self.collection = valid_edges.drop(drop_cols, axis=1)
        else:
            self.collection = self.collection.append(
                                valid_edges.drop(drop_cols,
                                                 axis=1), ignore_index=True)

        self.tell('current collection has', self.collection.shape[0], 'edges')
        # TODO: rename "query" to "protein" for consistency

    def write_graph(self, filename):
        self.collection.to_csv(filename, sep='\t')

    def compute_graph(self):
        collection_file = os.path.join(self.graphs_dir, self.alias)
        if not os.path.exists(collection_file):
            self.tell('Loading STRING core id')
            core_ids = []
            for line in open(self.core_ids, 'r'):
                if line != '':
                    core_ids.append(line.strip())

            # check whether the database fot the fasta file is created
            if not (os.path.exists(self.fasta + '.phr') and
                    os.path.exists(self.fasta + '.pin')):
                self.tell('Creating a blast database from the input fasta')
                fasta = self.fasta
                subprocess.call(f"makeblastdb -in {fasta} -out {fasta} " +
                                "-dbtype prot",
                                shell=True)

            self.tell('Computing Orthologs using', self.cpus, 'cores')
            params = []
            for core_id in core_ids:
                fasta = os.path.join(self.string_dir, core_id + '.faa')
                out = self.alias + '_AND_' + core_id
                params.append([self.fasta, self.db, fasta, fasta, out,
                               self.orthologs_dir, self.protein_format])
            with Pool(self.cpus) as p:
                p.starmap(compute_ortholog, params)

            # once the orthologs are computed, we transfer links from STRING
            should_process = {}
            graph_index = {}
            graph = {}
            first_line = True
            current_organism = -1
            for line in open(self.string_links):
                fields = line.strip().split()
                if first_line:
                    first_line = False

                    # identify the index of the interesting
                    # graphs in the STRING file
                    for i, field in enumerate(fields):
                        if field in self.interesting_graphs:
                            graph_index[field] = i
                else:
                    org_id = fields[0].split('.')[0]

                    # check if the current organism should be processed
                    if org_id not in should_process.keys():
                        should_process[org_id] =\
                            self.should_be_processed(org_id)
                        if not should_process[org_id]:
                            self.tell('Ignoring organism', org_id)
                    if not should_process[org_id]:
                        continue

                    if current_organism == -1:
                        current_organism = org_id
                        graph = self.clean_graph()
                    elif current_organism != org_id:
                        # make a pandas from the string graph
                        graph_df = pd.DataFrame.from_dict(graph)
                        # keeping only the relevant links
                        # (those with some information in at
                        # least one of the models)
                        condition = None
                        for g in self.interesting_graphs:
                            if condition is None:
                                condition = graph_df[g] > 0
                            else:
                                condition |= graph_df[g] > 0

                        # transfer links
                        self.process_graph(current_organism,
                                           graph_df[condition])

                        # clean graph and update current organism
                        graph = self.clean_graph()
                        current_organism = org_id

                    if should_process[org_id]:
                        graph['protein 1'].append(fields[0])
                        graph['protein 2'].append(fields[1])
                        for g in self.interesting_graphs:
                            graph[g].append(int(fields[graph_index[g]]))

            # we make sure we don't miss possible links from
            # the last organism in STRING
            if should_process[current_organism]:
                # make a pandas from the string graph
                graph_df = pd.DataFrame.from_dict(graph)
                # keeping only the relevant links (those with some
                # information in at least one of the models)
                condition = None
                for g in self.interesting_graphs:
                    if condition is None:
                        condition = graph_df[g] > 0
                    else:
                        condition |= graph_df[g] > 0

                # transfer links
                self.process_graph(current_organism, graph_df[condition])

            Graph.assert_lexicographical_order(self.collection,
                                               p1='query1', p2='query2')
            self.collection.to_pickle(os.path.join(collection_file))
        else:
            self.tell('Graph collection file found, skipping computation...')
            self.collection = pd.read_pickle(collection_file)
            Graph.assert_lexicographical_order(self.collection, p1='query1',
                                               p2='query2')
