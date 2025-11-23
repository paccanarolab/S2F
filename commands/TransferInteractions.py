import os
import pandas as pd
from Utils import ColourClass, Configuration, FancyApp, Utilities
from graphs.collection import Collection


class TransferInteractions(FancyApp.FancyApp):

    def __init__(self, args):
        super(TransferInteractions, self).__init__()
        self.run_config = os.path.expanduser(args.run_config)
        Configuration.load_run(self.run_config)
        run_conf = Configuration.RUN_CONFIG
        self.config_file = os.path.expanduser(
            run_conf.get('configuration', 'config_file'))
        Configuration.load_configuration(self.config_file)
        self.colour = ColourClass.bcolors.BOLD_GREEN
        self.alias = run_conf.get('configuration', 'alias')
        self.fasta = os.path.expanduser(run_conf.get('configuration',
                                                     'fasta'))
        self.installation_directory = os.path.expanduser(
            Configuration.CONFIG.get('directories', 'installation_directory'))
        self.string_links = Configuration.CONFIG.get('databases',
                                                     'string_links')
        self.string_dir = os.path.join(self.installation_directory,
                                       'data/STRINGSequences')
        self.core_ids = os.path.join(
            self.installation_directory, 'data/coreIds')
        self.orthologs_dir = os.path.join(
            self.installation_directory, 'orthologs')
        self.graphs_dir = os.path.join(
            self.installation_directory, 'graphs/collection')
        self.output_dir = os.path.join(self.installation_directory,
                                       'output', self.alias)
        self.max_evalue = args.max_evalue
        self.perc = args.perc
        self.positives = args.positives
        self.protein_format = run_conf.get('functions',
                                           'fasta_id_parser',
                                           fallback='uniprot')
        self.collection_chunk_size = run_conf.getint(
            'graphs', 'collection_chunk_size', fallback=200000)
        self.recompute_orthologs = run_conf.getboolean(
            'graphs', 'recompute_orthologs', fallback=True)
        self.cpu = os.path.expanduser(run_conf.get('configuration',
                                                   'cpu',
                                                   fallback='infer'))
        if self.cpu == 'infer':
            # https://docs.python.org/3/library/os.html#os.cpu_count
            self.cpu = len(os.sched_getaffinity(0))
        else:
            self.cpu = int(self.cpu)

        self.tell('Extracting list of proteins from fasta file')
        if self.protein_format == 'uniprot':
            fasta_id_parser = Utilities.keep_uniprot_accession
        else:
            fasta_id_parser = Utilities.keep_entire_prot_id
        self.proteins = Utilities.extract_indices_from_fasta(
            self.fasta,
            processing_func=fasta_id_parser
        )
        self.proteins.to_pickle(os.path.join(self.output_dir, 'proteins.df'))

    def run(self):
        col = Collection(self.fasta, self.proteins, self.string_dir,
                         self.string_links, self.core_ids, self.output_dir,
                         self.orthologs_dir, self.graphs_dir, self.alias,
                         self.cpu, None, self.max_evalue, self.perc,
                         self.positives, self.protein_format,
                         recompute_orthologs=self.recompute_orthologs,
                         chunk_size=self.collection_chunk_size
                         )
        col.compute_graph()
