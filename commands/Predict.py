from Utils import *
from seeds import interpro
from seeds import hmmer
from graphs import Graph
from graphs import collection
from graphs import homology
from graphs import combination
from GOTool import GeneOntology

from scipy import sparse
from scipy.sparse.linalg import lsqr
from sklearn.preprocessing import binarize

import pandas as pd
import numpy as np
import os
import subprocess
import sys


class Predict(FancyApp.FancyApp):

    def __init__(self, args):
        super(Predict, self).__init__()
        self.colour = ColourClass.bcolors.OKGREEN
        self.run_config = os.path.expanduser(args.run_config)
        # load configurations
        if os.path.exists(self.run_config):
            self.tell('Run configuration provided, loading: ', self.run_config)
            Configuration.load_run(self.run_config)

            self.config_file = os.path.expanduser(Configuration.RUN_CONFIG.get('configuration', 'config_file'))
            self.alias = Configuration.RUN_CONFIG.get('configuration', 'alias')
            self.obo = os.path.expanduser(Configuration.RUN_CONFIG.get('configuration', 'obo'))
            self.fasta = os.path.expanduser(Configuration.RUN_CONFIG.get('configuration', 'fasta'))
            self.cpu = os.path.expanduser(Configuration.RUN_CONFIG.get('configuration', 'cpu', fallback='infer'))

            self.combined_graph = Configuration.RUN_CONFIG.get('graphs', 'combined_graph', fallback='compute')
            self.graph_collection = Configuration.RUN_CONFIG.get('graphs', 'graph_collection', fallback='compute')
            self.homology_graph = Configuration.RUN_CONFIG.get('graphs', 'homology_graph', fallback='compute')

            self.interpro_output = Configuration.RUN_CONFIG.get('seeds', 'interpro_output', fallback='compute')
            self.hmmer_output = Configuration.RUN_CONFIG.get('seeds', 'hmmer_output', fallback='compute')

            self.hmmer_blacklist = Configuration.RUN_CONFIG.get('blacklists', 'hmmer_blacklist', fallback='compute')
            self.transfer_blacklist = Configuration.RUN_CONFIG.get('blacklists',
                                                                   'transfer_blacklist', fallback='compute')

            self.goa_clamp = Configuration.RUN_CONFIG.get('functions', 'goa_clamp', fallback='compute')

        else:
            Configuration.load_configuration(self.config_file)

            self.config_file = os.path.expanduser(args.config_file)
            self.alias = args.alias
            self.obo = os.path.expanduser(args.obo)
            self.fasta = os.path.expanduser(args.fasta)
            self.cpu = os.path.expanduser(args.cpu)

            self.combined_graph = os.path.expanduser(args.combined_graph)
            self.graph_collection = os.path.expanduser(args.graph_collection)
            self.homology_graph = os.path.expanduser(args.homology_graph)

            self.interpro_output = os.path.expanduser(args.interpro_output)
            self.hmmer_output = os.path.expanduser(args.hmmer_seed)

            self.hmmer_blacklist = os.path.expanduser(args.hmmer_blacklist)
            self.transfer_blacklist = os.path.expanduser(args.transfer_blacklist)

            self.goa_clamp = os.path.expanduser(args.goa_clamp)

        self.installation_directory = os.path.expanduser(
            Configuration.CONFIG.get('directories', 'installation_directory'))

        self.interpro = Configuration.CONFIG.get('commands', 'interpro')
        self.hmmer = Configuration.CONFIG.get('commands', 'hmmer')
        self.blastp = Configuration.CONFIG.get('commands', 'blastp')
        self.makeblastdb = Configuration.CONFIG.get('commands', 'makeblastdb')

        self.string_links = Configuration.CONFIG.get('databases', 'string_links')
        self.string_sequences = Configuration.CONFIG.get('databases', 'string_sequences')
        self.string_species = Configuration.CONFIG.get('databases', 'string_species')
        self.uniprot_sprot = Configuration.CONFIG.get('databases', 'uniprot_sprot')
        self.uniprot_goa = Configuration.CONFIG.get('databases', 'uniprot_goa')
        self.filtered_goa = Configuration.CONFIG.get('databases', 'filtered_goa')
        self.filtered_sprot = Configuration.CONFIG.get('databases', 'filtered_sprot')

        ec = Configuration.CONFIG.get('options', 'evidence_codes')
        self.evidence_codes = ec if ec == 'experimental' else ec.split(',')

        if (not os.path.exists(self.config_file) or
                self.alias == '' or
                not os.path.exists(self.obo) or
                not os.path.exists(self.fasta)):
            self.warning('The configuration file must include readable files and an alias, please revise these'
                         'arguments:\n'
                         '\tconfig_file\n'
                         '\talias\n'
                         '\tobo\n'
                         '\tfasta')
            sys.exit(1)

        self.summary_and_continue()

        self.output_dir = os.path.join(self.installation_directory, 'output', self.alias)
        self.seed_dir_IP = os.path.join(self.installation_directory, 'seeds/interpro')
        self.seed_dir_H = os.path.join(self.installation_directory, 'seeds/hmmer')
        self.go = GeneOntology.GeneOntology(self.obo, verbose=True)
        self.terms = None
        self.proteins = None

        if self.cpu == 'infer':
            self.cpu = len(os.sched_getaffinity(0))  # https://docs.python.org/3/library/os.html#os.cpu_count

    def run(self):

        # 1. create output directory
        self.create_output_directory()
        # 2. Exract a list of terms and proteins
        self.make_indices()
        # 3. run (or reuse) interpro
        ip_seed = self.run_interpro()
        # 4. run (or reuse) hmmer
        hmmer_seed = self.run_hmmer()
        # 5. run (or reuse) graph collection
        graph_collection = self.run_graph_collection()
        for k, g in graph_collection.items():
            sparse.save_npz(os.path.join(self.output_dir, k+'_collection.npz'), g)
        # 5. run (or reuse) homology graph
        graph_homology = self.run_graph_homology()
        sparse.save_npz(os.path.join(self.output_dir, 'homology.npz'), graph_homology)
        # 6. GOA clamp if provided
        # 7. graph combination
        combined_graph = self.combine_graphs(graph_collection, graph_homology, ip_seed)
        # 8. diffuse interpro
        # 9. diffuse hmmer
        # 10. output combination
        pass

    def create_output_directory(self):
        if not os.path.isdir(self.output_dir):
            self.tell('Creating output directory')
            os.mkdir(self.output_dir)
        self.tell(self.output_dir, 'will be used as output directory')

    def make_indices(self):
        self.tell('Extracting list of proteins from fasta file')
        self.proteins = Utilities.extract_indices_from_fasta(self.fasta)
        self.proteins.to_pickle(os.path.join(self.output_dir, 'proteins.df'))

        self.tell('Building GO structure')
        self.go.build_structure()
        self.tell('Extracting list of GO Terms from structure')
        self.terms = pd.DataFrame(list(enumerate(sorted(self.go.terms.keys()))))
        self.terms.columns = ['term idx', 'term id']
        self.terms.set_index('term id', inplace=True)
        self.terms.to_pickle(os.path.join(self.output_dir, 'terms.df'))

    def run_interpro(self):
        if not os.path.exists(self.interpro_output) or self.interpro_output == 'compute':
            self.interpro_output = os.path.join(self.output_dir, self.alias + '_IP')
            if not os.path.exists(self.interpro_output):
                self.tell('Running InterPro')
                command = self.interpro + ' -i ' + self.fasta + ' -goterms -iprlookup -f TSV -o ' + self.interpro_output
                subprocess.call(command, shell=True)
            else:
                self.tell('InterPro file found, skipping computation...')
        return self.combine_interpro()

    def combine_interpro(self):
        seed_file = os.path.join(self.seed_dir_IP, self.alias + '.seed.npz')
        if not os.path.exists(seed_file):
            self.tell('Building InterPro seed file')
            interpro_seed = interpro.InterProSeed(self.interpro_output, self.proteins, self.terms, self.go)
            # TODO: methods is here only to debug stuff, remove
            methods = interpro_seed.process_output()
            for k, s in methods.items():
                s.to_pickle(os.path.join(self.output_dir, 'ipseed_'+k))
            seed = interpro_seed.get_seed()
            sparse.save_npz(seed_file, seed)
        else:
            self.tell('InterPro seed file found')
            seed = sparse.load_npz(seed_file)
        return seed

    def run_hmmer(self):
        out = os.path.join(self.output_dir, self.alias + '.hmmer')
        if not os.path.exists(self.hmmer_output) or self.hmmer_output == 'compute':
            self.hmmer_output = os.path.join(self.output_dir, self.alias + '_H')
            if not os.path.exists(self.hmmer_output):
                self.tell('Running HMMer')
                command = self.hmmer + ' --cpu ' + str(self.cpu) + ' --noali -o ' + out + \
                    ' --tblout ' + self.hmmer_output + ' ' + self.fasta + ' ' + self.filtered_sprot
                self.tell(command)
                subprocess.call(command, shell=True)
            else:
                self.tell('HMMer file found, skipping computation...')
        return self.build_hmeer_seed()

    def build_hmeer_seed(self):
        seed_file = os.path.join(self.seed_dir_H, self.alias + '.seed.npz')
        evalue_file = os.path.join(self.seed_dir_H, self.alias + '.evalue')
        if not os.path.exists(seed_file):
            self.tell('Building HMMer seed file')
            hmmer_seed = hmmer.HMMerSeed(self.hmmer_output, self.proteins, self.terms, self.go,
                                         self.hmmer_blacklist, self.filtered_goa)
            hmmer_seed.process_output(evalue_file=evalue_file)
            seed = hmmer_seed.get_seed()
            sparse.save_npz(seed_file, seed)
        else:
            self.tell('HMMer seed file found, skipping computation...')
            seed = sparse.load_npz(seed_file)
        return seed

    def run_graph_collection(self):
        string_dir = os.path.join(self.installation_directory, 'data/STRINGSequences')
        core_ids = os.path.join(self.installation_directory, 'data/coreIds')
        orthologs_dir = os.path.join(self.installation_directory, 'orthologs')
        graphs_dir = os.path.join(self.installation_directory, 'graphs/collection')
        col = collection.Collection(self.fasta, self.proteins, string_dir, self.string_links, core_ids,
                                    self.output_dir, orthologs_dir, graphs_dir, self.alias, self.cpu,
                                    self.transfer_blacklist, 1e-6, 80.0, 60.0)
        col.compute_graph()
        return col.get_graph()

    def run_graph_homology(self):
        graphs_dir = os.path.join(self.installation_directory, 'graphs/homology')
        hom = homology.Homology(self.fasta, self.proteins, graphs_dir, self.alias)
        hom.compute_graph()
        return hom.get_graph()

    def goa_clamp(self):
        pass

    def combine_graphs(self, graph_collection, graph_homology, ip_seed):
        comb = combination.Combination(self.proteins, graph_collection, graph_homology, ip_seed)
        comb.compute_graph()
        return comb.get_graph()

    def diffuse(self):
        pass

    def combine_diffusion(self):
        pass

    def summary_and_continue(self):
        summary = 'These are the loaded values \n'
        summary += ColourClass.coloured_string(ColourClass.bcolors.HEADER, '\t[INSTALLATION CONFIGURATION]\n')
        summary += '\tInstallation directory:\t\t' + self.installation_directory + '\n\n'

        summary += '\tPath to InterPro executable:\t' + self.interpro + '\n'
        summary += '\tPath to HMMer:\t\t\t' + self.hmmer + '\n'
        summary += '\tPath to blastp:\t\t\t' + self.blastp + '\n'
        summary += '\tPath to makeblastdb:\t\t' + self.makeblastdb + '\n\n'

        # TODO: These warnings are extremely unlikely to occur because the installation should put all of these in the
        # right values. This code could be simplified.
        dbs = [
            ('STRING interaction database:\t', self.string_links),
            ('STRING sequences database:\t', self.string_sequences),
            ('STRING species database:\t', self.string_species),
            ('UniProt SwissProt:\t\t', self.uniprot_sprot),
            ('UniProt GOA:\t\t\t', self.uniprot_goa)
        ]
        for desc, database in dbs:
            summary += '\t' + desc
            if database == 'download':
                summary += FancyApp.FancyApp.warning_text(database) + '\n'
            else:
                summary += database + '\n'

        summary += '\n\tFiltered UniProt GOA:\t\t' + self.filtered_goa + '\n'
        summary += '\tFiltered UniProt SwissProt:\t' + self.filtered_sprot + '\n\n'

        summary += '\tEvidence Codes:\t\t\t' + str(self.evidence_codes) + '\n\n'

        summary += ColourClass.coloured_string(ColourClass.bcolors.HEADER, '\t[RUN CONFIGURATION]\n')
        # summary += '\tInstallation configuration file: ' + self.config_file + '\n'
        summary += '\tAlias:\t\t\t\t' + self.alias + '\n'
        summary += '\tOBO file:\t\t\t' + self.obo + '\n'
        summary += '\tFasta file:\t\t\t' + self.fasta + '\n\n'

        summary += '\tCombined Graph:\t\t\t' + self.combined_graph + '\n'
        summary += '\tGraph Collection:\t\t' + self.graph_collection + '\n'
        summary += '\tHomology Graph:\t\t\t' + self.homology_graph + '\n\n'

        summary += '\tInterPro seed:\t\t\t' + self.interpro_output + '\n'
        summary += '\tHMMer seed:\t\t\t' + self.hmmer_output + '\n'

        summary += '\tHMMer black list:\t\t' + self.hmmer_blacklist + '\n'
        summary += '\tTransfer black list:\t\t' + self.transfer_blacklist + '\n'

        summary += '\tGOA clamp:\t\t\t' + self.goa_clamp + '\n'

        summary += 'Do you want to continue with these settings?'

        if Utilities.query_yes_no(summary):
            return
        sys.exit()
