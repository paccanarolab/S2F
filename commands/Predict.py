import os
import subprocess
import sys

import pandas as pd
import numpy as np
from scipy import sparse

from diffusion import Diffusion
from diffusion.S2FLabelPropagation import S2FLabelPropagation
from GOTool import GeneOntology
from graphs import collection, combination, homology, Graph
from seeds import hmmer, interpro
from Utils import ColourClass, Configuration, FancyApp, Utilities


class Predict(FancyApp.FancyApp):

    def __init__(self, args):
        super(Predict, self).__init__()
        self.colour = ColourClass.bcolors.OKGREEN
        if args.run_config != 'arguments':
            self.run_config = os.path.expanduser(args.run_config)
        else:
            self.run_config = 'arguments'
        # load configurations
        if os.path.exists(self.run_config):
            self.tell('Run configuration provided, loading: ', self.run_config)
            Configuration.load_run(self.run_config)
            run_conf = Configuration.RUN_CONFIG

            self.config_file = os.path.expanduser(
                run_conf.get('configuration', 'config_file'))
            self.alias = run_conf.get('configuration', 'alias')
            self.obo = os.path.expanduser(run_conf.get('configuration',
                                                       'obo'))
            self.fasta = os.path.expanduser(run_conf.get('configuration',
                                                         'fasta'))
            self.cpu = os.path.expanduser(run_conf.get('configuration',
                                                       'cpu',
                                                       fallback='infer'))

            self.combined_graph = run_conf.get('graphs', 'combined_graph',
                                               fallback='compute')
            self.graph_collection = run_conf.get('graphs', 'graph_collection',
                                                 fallback='compute')
            self.homology_graph = run_conf.get('graphs', 'homology_graph',
                                               fallback='compute')

            self.interpro_output = run_conf.get('seeds',
                                                'interpro_output',
                                                fallback='compute')
            self.hmmer_output = run_conf.get('seeds',
                                             'hmmer_output',
                                             fallback='compute')

            self.hmmer_blacklist = run_conf.get('blacklists',
                                                'hmmer_blacklist',
                                                fallback='compute')
            self.transfer_blacklist = run_conf.get('blacklists',
                                                   'transfer_blacklist',
                                                   fallback='compute')

            self.goa_clamp = run_conf.get('functions',
                                          'goa_clamp',
                                          fallback='compute')
            self.unattended_mode = run_conf.getboolean('functions',
                                                       'unattended',
                                                       fallback=False)
            self.write_collection = run_conf.getboolean('functions',
                                                        'write_collection',
                                                        fallback=False)
            self.fasta_id_parser = run_conf.get('functions',
                                                'fasta_id_parser',
                                                fallback='uniprot')
        else:
            self.config_file = os.path.expanduser(args.config_file)

            Configuration.load_configuration(self.config_file)

            self.alias = args.alias
            self.obo = os.path.expanduser(args.obo)
            self.fasta = os.path.expanduser(args.fasta)
            self.cpu = os.path.expanduser(args.cpu)

            self.combined_graph = os.path.expanduser(args.combined_graph)
            self.graph_collection = os.path.expanduser(args.graph_collection)
            self.homology_graph = os.path.expanduser(args.homology_graph)

            self.interpro_output = os.path.expanduser(args.interpro_output)
            self.hmmer_output = os.path.expanduser(args.hmmer_output)

            self.hmmer_blacklist = os.path.expanduser(args.hmmer_blacklist)
            self.transfer_blacklist = os.path.expanduser(
                                            args.transfer_blacklist)

            self.goa_clamp = os.path.expanduser(args.goa_clamp)
            self.unattended_mode = args.unattended is True
            self.write_collection = args.write_collection is True
            self.fasta_id_parser = args.fasta_id_parser

        self.installation_directory = os.path.expanduser(
            Configuration.CONFIG.get('directories', 'installation_directory'))

        self.interpro = Configuration.CONFIG.get('commands', 'interpro')
        self.hmmer = Configuration.CONFIG.get('commands', 'hmmer')
        self.blastp = Configuration.CONFIG.get('commands', 'blastp')
        self.makeblastdb = Configuration.CONFIG.get('commands', 'makeblastdb')

        self.string_links = Configuration.CONFIG.get('databases',
                                                     'string_links')
        self.string_sequences = Configuration.CONFIG.get('databases',
                                                         'string_sequences')
        self.string_species = Configuration.CONFIG.get('databases',
                                                       'string_species')
        self.string_core_only = Configuration.CONFIG.get('databases',
                                                         'string_core_only')
        self.uniprot_sprot = Configuration.CONFIG.get('databases',
                                                      'uniprot_sprot')
        self.uniprot_goa = Configuration.CONFIG.get('databases',
                                                    'uniprot_goa')
        self.filtered_goa = Configuration.CONFIG.get('databases',
                                                     'filtered_goa')
        self.filtered_sprot = Configuration.CONFIG.get('databases',
                                                       'filtered_sprot')

        ec = Configuration.CONFIG.get('options', 'evidence_codes')
        self.evidence_codes = ec if ec == 'experimental' else ec.split(',')

        if (not os.path.exists(self.config_file) or
                self.alias == '' or
                not os.path.exists(self.obo) or
                not os.path.exists(self.fasta)):
            self.warning("The configuration file must include readable " +
                         "files and an alias, please revise these" +
                         'arguments:\n' +
                         '\tconfig_file\n' +
                         '\talias\n' +
                         '\tobo\n' +
                         '\tfasta')
            sys.exit(1)

        if not self.unattended_mode:
            self.summary_and_continue()

        self.output_dir = os.path.join(self.installation_directory,
                                       'output', self.alias)
        self.seed_dir_IP = os.path.join(self.installation_directory,
                                        'seeds/interpro')
        self.seed_dir_H = os.path.join(self.installation_directory,
                                       'seeds/hmmer')
        self.combination_dir = os.path.join(self.installation_directory,
                                            'graphs', 'combined')
        self.go = GeneOntology.GeneOntology(self.obo, verbose=True)
        self.terms = None
        self.proteins = None
        self.prediction = None
        self.clamp_matrix = None

        if self.cpu == 'infer':
            # https://docs.python.org/3/library/os.html#os.cpu_count
            self.cpu = len(os.sched_getaffinity(0))
        else:
            self.cpu = int(self.cpu)

    def run(self):
        # 1. create output directory
        self.create_output_directory()
        # 2. Exract a list of terms and proteins
        self.make_indices()
        # 3. run (or reuse) interpro
        self.check_progress()
        if self.load_ip_seed:
            ip_seed = self.run_interpro()
        else:
            ip_seed = None
        # 4. run (or reuse) hmmer
        if self.load_hmmer_seed:
            hmmer_seed = self.run_hmmer()
        else:
            hmmer_seed = None
        # 5. run (or reuse) graph collection
        if self.load_graph_collection:
            graph_collection = self.run_graph_collection()
            if self.write_collection:
                for k, g in graph_collection.items():
                    sparse.save_npz(os.path.join(self.output_dir,
                                                 k+'_collection.npz'),
                                    g)
        else:
            graph_collection = None
        # 6. run (or reuse) homology graph
        if self.load_homology_graph:
            graph_homology = self.run_graph_homology()
            sparse.save_npz(os.path.join(self.output_dir, 'homology.npz'),
                            graph_homology)
        else:
            graph_homology = None
        # 7. GOA clamp if provided
        if os.path.exists(self.goa_clamp):
            self.load_clamp()
            self.tell('Clamping InterPro seed...')
            ip_seed = self.clamp(ip_seed)
            self.tell('Clamping HMMER seed...')
            hmmer_seed = self.clamp(hmmer_seed)
        # 8. graph combination
        if self.load_combined_graph:
            combined_graph = self.combine_graphs(graph_collection,
                                                 graph_homology,
                                                 ip_seed)
        else:
            combined_graph = None
        # 9. prepare diffusion kernel
        if self.calculate_diffusion:
            kernel_params = {'lambda': 0.1}
            diff = S2FLabelPropagation(combined_graph,
                                       self.proteins,
                                       self.terms)
            diff.compute_kernel(**kernel_params)
            # 10. diffuse interpro
            self.tell('Diffusing InterPro seed')
            ip_diff = diff.diffuse(ip_seed)
            sparse.save_npz(self.ip_diff_file + '.npz', ip_diff)
            diff.write_results(self.ip_diff_file)
            # 11. diffuse hmmer
            self.tell('Diffusing HMMER seed')
            hmmer_diff = diff.diffuse(hmmer_seed)
            sparse.save_npz(self.hmmer_diff_file + '.npz', hmmer_diff)
            diff.write_results(self.hmmer_diff_file)
        else:
            self.tell('Diffusion files found, skipping computation...')
            ip_diff = sparse.load_npz(self.ip_diff_file + '.npz')
            hmmer_diff = sparse.load_npz(self.hmmer_diff_file + '.npz')
        # 12. output combination
        if os.path.exists(self.goa_clamp):
            self.tell('Clamping diffused InterPro')
            ip_diff = self.clamp(ip_diff)
            self.tell('Clamping diffused HMMER')
            hmmer_diff = self.clamp(hmmer_diff)
        self.tell('Adding HMMER and InterPro diffusions')
        self.prediction = 0.9 * ip_diff + 0.1 * hmmer_diff
        self.prediction = sparse.coo_matrix(self.prediction)
        self.tell('Saving prediction to file')
        self.write_prediction(os.path.join(self.output_dir, 'prediction.df'))

    def check_progress(self):
        self.ip_seed_file = os.path.join(
            self.seed_dir_IP, self.alias + '.seed.npz')
        self.ip_diff_file = os.path.join(
            self.output_dir, 'ip_seed.diffusion')

        self.hmmer_seed_file = os.path.join(
            self.seed_dir_H, self.alias + '.seed.npz')
        self.hmmer_evalue_file = os.path.join(
            self.seed_dir_H, self.alias + '.evalue')
        self.hmmer_diff_file = os.path.join(
            self.output_dir, 'hmmer_seed.diffusion')

        self.string_dir = os.path.join(self.installation_directory,
                                       'data/STRINGSequences')
        self.core_ids = os.path.join(
            self.installation_directory, 'data/coreIds')
        self.orthologs_dir = os.path.join(
            self.installation_directory, 'orthologs')
        self.graphs_dir = os.path.join(
            self.installation_directory, 'graphs/collection')

        self.combined_graph_filename = os.path.join(self.combination_dir,
                                                    self.alias)
        self.combined_graph_sparse_filename =\
            f"{self.combined_graph_filename}.npz"

        # we assume everything was done
        self.load_ip_seed = False
        self.load_hmmer_seed = False
        self.load_graph_collection = False
        self.load_homology_graph = False
        self.load_combined_graph = False
        self.calculate_diffusion = False
        # the strategy is to check from latest to earliest in order to load
        # only what we need.
        if os.path.exists(self.ip_diff_file) and os.path.exists(
                self.hmmer_diff_file):
            self.calculate_diffusion = False
            return
        else:
            self.calculate_diffusion = True
            self.load_ip_seed = True
            self.load_hmmer_seed = True
            self.load_combined_graph = True
            if os.path.exists(self.combined_graph_filename):
                return
            else:
                self.load_homology_graph = True
                self.load_graph_collection = True

    def write_prediction(self, filename):
        Diffusion._write_results(self.prediction, self.proteins, self.terms,
                                 filename)

    def create_output_directory(self):
        if not os.path.isdir(self.output_dir):
            self.tell('Creating output directory')
            os.mkdir(self.output_dir)
        self.tell(self.output_dir, 'will be used as output directory')

    def make_indices(self):
        self.tell('Extracting list of proteins from fasta file')
        if self.fasta_id_parser == 'uniprot':
            fasta_id_parser = Utilities.keep_uniprot_accession
        else:
            fasta_id_parser = Utilities.keep_entire_prot_id
        self.proteins = Utilities.extract_indices_from_fasta(
            self.fasta,
            processing_func=fasta_id_parser
        )
        self.proteins.to_pickle(os.path.join(self.output_dir, 'proteins.df'))

        self.tell('Building GO structure')
        self.go.build_structure()
        self.tell('Extracting list of GO Terms from structure')
        self.terms = pd.DataFrame(
                        list(enumerate(sorted(self.go.terms.keys()))))
        self.terms.columns = ['term idx', 'term id']
        self.terms.set_index('term id', inplace=True)
        self.terms.to_pickle(os.path.join(self.output_dir, 'terms.df'))

    def run_interpro(self):
        if not os.path.exists(self.interpro_output) or\
                self.interpro_output == 'compute':
            self.interpro_output = os.path.join(self.output_dir,
                                                self.alias + '_IP')
            if not os.path.exists(self.interpro_output):
                self.tell('Running InterPro')
                command = self.interpro + ' -i ' + self.fasta +\
                                          ' -goterms -iprlookup -f TSV -o ' +\
                                          self.interpro_output
                subprocess.call(command, shell=True)
            else:
                self.tell('InterPro file found, skipping computation...')
        return self.combine_interpro()

    def combine_interpro(self):
        if not os.path.exists(self.ip_seed_file):
            self.tell('Building InterPro seed file')
            interpro_seed = interpro.InterProSeed(self.interpro_output,
                                                  self.proteins,
                                                  self.terms, self.go,
                                                  self.fasta_id_parser
                                                  )
            # TODO: methods is here only to debug stuff, remove
            methods = interpro_seed.process_output()
            for k, s in methods.items():
                s.to_pickle(os.path.join(self.output_dir, 'ipseed_'+k))
            seed = interpro_seed.get_seed()
            sparse.save_npz(self.ip_seed_file, seed)
        else:
            self.tell('InterPro seed file found')
            seed = sparse.load_npz(self.ip_seed_file)
        return seed

    def run_hmmer(self):
        out = os.path.join(self.output_dir, self.alias + '.hmmer')
        if not os.path.exists(self.hmmer_output) or\
                self.hmmer_output == 'compute':
            self.hmmer_output = os.path.join(self.output_dir,
                                             self.alias + '_H')
            if not os.path.exists(self.hmmer_output):
                self.tell('Running HMMer')
                command = self.hmmer + ' --cpu ' + str(self.cpu) +\
                                       ' --noali -o ' + out + \
                                       ' --tblout ' + self.hmmer_output +\
                                       ' ' + self.fasta + ' ' +\
                                       self.filtered_sprot
                self.tell(command)
                subprocess.call(command, shell=True)
            else:
                self.tell('HMMer file found, skipping computation...')
        return self.build_hmeer_seed()

    def build_hmeer_seed(self):
        if not os.path.exists(self.hmmer_seed_file):
            self.tell('Building HMMer seed file')
            hmmer_seed = hmmer.HMMerSeed(self.hmmer_output, self.proteins,
                                         self.terms, self.go,
                                         self.get_hmmer_blacklist(),
                                         self.filtered_goa, 
                                         self.fasta_id_parser
                                         )
            hmmer_seed.process_output(evalue_file=self.hmmer_evalue_file)
            seed = hmmer_seed.get_seed()
            sparse.save_npz(self.hmmer_seed_file, seed)
        else:
            self.tell('HMMer seed file found, skipping computation...')
            seed = sparse.load_npz(self.hmmer_seed_file)
        return seed

    def get_hmmer_blacklist(self):
        if self.hmmer_blacklist != 'compute':
            return pd.read_csv(
                self.hmmer_blacklist,
                header=None,
                names=['tax_id']
            ).tax_id.astype('str').tolist()
        return None

    def run_graph_collection(self):
        col = collection.Collection(self.fasta, self.proteins, self.string_dir,
                                    self.string_links, self.core_ids,
                                    self.output_dir, self.orthologs_dir,
                                    self.graphs_dir, self.alias, self.cpu,
                                    self.get_transfer_blacklist(),
                                    1e-6, 80.0, 60.0, self.fasta_id_parser)
        col.compute_graph()
        return col.get_graph()

    def get_transfer_blacklist(self):
        if self.hmmer_blacklist != 'compute':
            return pd.read_csv(
                self.transfer_blacklist,
                header=None,
                names=['tax_id']
            ).tax_id.astype('str').tolist()
        return None

    def run_graph_homology(self):
        graphs_dir = os.path.join(self.installation_directory,
                                  'graphs/homology')
        hom = homology.Homology(self.fasta, self.proteins,
                                graphs_dir,
                                self.alias, self.fasta_id_parser)
        hom.compute_graph()
        return hom.get_graph()

    def load_clamp(self):
        self.tell('Loading GOA annotations to clamp...')
        self.go.load_clamp_file(self.goa_clamp, 'clamp')
        self.tell('Up-propagating clamp annotations...')
        self.go.up_propagate_annotations('clamp')
        self.tell('Building clamp matrix...')
        clamp_df = self.go.get_annotations('clamp')
        clamp_df = clamp_df.merge(self.proteins,
                                  left_on='Protein',
                                  right_index=True)
        clamp_df = clamp_df.merge(self.terms,
                                  left_on='GO ID',
                                  right_index=True)

        p_idx = clamp_df['protein idx'].values
        go_idx = clamp_df['term idx'].values

        self.clamp_matrix = sparse.coo_matrix((np.ones(len(clamp_df)),
                                               (p_idx, go_idx)),
                                              shape=(len(self.proteins),
                                                     len(self.terms)))

    def clamp(self, seeds):
        clamp_bigger = self.clamp_matrix > seeds
        return (seeds - seeds.multiply(clamp_bigger) +
                self.clamp_matrix.multiply(clamp_bigger))

    def combine_graphs(self, graph_collection, graph_homology, ip_seed):
        if not os.path.exists(self.combined_graph_filename):
            comb = combination.Combination(self.proteins, graph_collection,
                                           graph_homology, ip_seed)
            comb.compute_graph()
            comb.write_graph(self.combined_graph_filename)
            graph = comb.get_graph()
            sparse.save_npz(self.combined_graph_sparse_filename, graph)
        else:
            self.tell('found combined graph, skipping comptuation')
            graph = sparse.load_npz(self.combined_graph_sparse_filename)
        graph = Graph.fill_lower_triangle(graph)
        return graph

    def summary_and_continue(self):
        summary = 'These are the loaded values \n'
        summary += ColourClass.coloured_string(
                        ColourClass.bcolors.HEADER,
                        "\t[INSTALLATION CONFIGURATION]\n")
        summary += '\tInstallation directory:\t\t'
        summary += self.installation_directory + '\n\n'
        summary += '\tPath to InterPro executable:\t' + self.interpro + '\n'
        summary += '\tPath to HMMer:\t\t\t' + self.hmmer + '\n'
        summary += '\tPath to blastp:\t\t\t' + self.blastp + '\n'
        summary += '\tPath to makeblastdb:\t\t' + self.makeblastdb + '\n\n'

        # TODO: These warnings are extremely unlikely to occur because
        # the installation should put all of these in the
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
        summary += '\tFiltered UniProt SwissProt:\t' + self.filtered_sprot
        summary += '\n\n'

        summary += '\tEvidence Codes:\t\t\t' + str(self.evidence_codes)
        summary += '\n\n'

        summary += ColourClass.coloured_string(ColourClass.bcolors.HEADER,
                                               '\t[RUN CONFIGURATION]\n')
        # summary += '\tInstallation configuration file: ' +\
        #            self.config_file + '\n'
        summary += '\tAlias:\t\t\t\t' + self.alias + '\n'
        summary += '\tOBO file:\t\t\t' + self.obo + '\n'
        summary += '\tFasta file:\t\t\t' + self.fasta + '\n\n'

        summary += '\tCombined Graph:\t\t\t' + self.combined_graph + '\n'
        summary += '\tGraph Collection:\t\t' + self.graph_collection + '\n'
        summary += '\tHomology Graph:\t\t\t' + self.homology_graph + '\n\n'

        summary += '\tInterPro seed:\t\t\t' + self.interpro_output + '\n'
        summary += '\tHMMer seed:\t\t\t' + self.hmmer_output + '\n'

        summary += '\tHMMer black list:\t\t' + self.hmmer_blacklist + '\n'
        summary += '\tTransfer black list:\t\t'
        summary += self.transfer_blacklist + '\n'
        summary += '\tGOA clamp:\t\t\t' + self.goa_clamp + '\n'
        summary += '\tFASTA handling:\t\t\t' + self.fasta_id_parser + '\n'

        summary += 'Do you want to continue with these settings?'

        if Utilities.query_yes_no(summary):
            return
        sys.exit()
