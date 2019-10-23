from GOTool import GeneOntology
from Utils import *
from requests_html import HTMLSession

import os
import subprocess
import re
import sys


class Install(FancyApp.FancyApp):

    def __init__(self, args):
        super(Install, self).__init__()
        self.colour = ColourClass.bcolors.OKGREEN
        self.config_file = os.path.expanduser(args.config_file)
        # load configuration
        if os.path.exists(self.config_file):
            self.tell('Configuration file found, loading:', self.config_file)
            Configuration.load_configuration(self.config_file)

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
        else:
            self.tell('Creating configuration file: ', self.config_file)
            self.tell('Provided arguments will override default ones')

            # collect parameters from provided arguments
            self.installation_directory = os.path.expanduser(args.installation_directory)

            self.interpro = args.interpro
            self.hmmer = args.hmmer
            self.blastp = args.blastp
            self.makeblastdb = args.makeblastdb

            self.string_links = args.string_links
            self.string_sequences = args.string_sequences
            self.string_species = args.string_species
            self.uniprot_sprot = args.uniprot_swissprot
            self.uniprot_goa = args.uniprot_goa
            # These parameters have to be inferred in the installation
            self.filtered_goa = 'infer'
            self.filtered_sprot = 'infer'

            self.evidence_codes = args.evidence_codes

            Configuration.CONFIG.add_section('directories')
            Configuration.CONFIG.set('directories', 'installation_directory', self.installation_directory)

            Configuration.CONFIG.add_section('commands')
            Configuration.CONFIG.set('commands', 'interpro', self.interpro)
            Configuration.CONFIG.set('commands', 'hmmer', self.hmmer)
            Configuration.CONFIG.set('commands', 'blastp', self.blastp)
            Configuration.CONFIG.set('commands', 'makeblastdb', self.makeblastdb)

            Configuration.CONFIG.add_section('databases')
            Configuration.CONFIG.set('databases', 'string_links', self.string_links)
            Configuration.CONFIG.set('databases', 'string_sequences', self.string_sequences)
            Configuration.CONFIG.set('databases', 'string_species', self.string_species)
            Configuration.CONFIG.set('databases', 'uniprot_sprot', self.uniprot_sprot)
            Configuration.CONFIG.set('databases', 'uniprot_goa', self.uniprot_goa)

            Configuration.CONFIG.add_section('options')
            if self.evidence_codes == 'experimental':
                self.evidence_codes = GeneOntology.GeneOntology.EXPERIMENTAL_EVIDENCE_CODES
            Configuration.CONFIG.set('options', 'evidence_codes', ','.join(self.evidence_codes))
        self.summary_and_continue()

    def run(self):
        self.create_directories()
        self.download_databases()
        self.process_files()
        self.tell('Updating configuration file')
        with open(self.config_file, 'w') as cf:
            Configuration.CONFIG.write(cf)

    def create_directories(self):
        self.tell('Creating directories')
        self.create_directory('')
        self.create_directory('data')
        self.create_directory('data/STRINGSequences')
        self.create_directory('data/UniprotKB')
        self.create_directory('data/UserSequences')
        self.create_directory('graphs')
        self.create_directory('graphs/collection')
        self.create_directory('graphs/combined')
        self.create_directory('graphs/homology')
        self.create_directory('INPUT')
        self.create_directory('INPUT/fastas')
        self.create_directory('INPUT/targetFiles')
        self.create_directory('orthologs')
        self.create_directory('output')
        self.create_directory('seeds')
        self.create_directory('seeds/hmmer')
        self.create_directory('seeds/interpro')

    def create_directory(self, directory):
        target = os.path.expanduser(os.path.join(self.installation_directory, directory))
        if not os.path.exists(target):
            os.mkdir(target)

    def download_databases(self):
        # STRING
        s = HTMLSession()
        r = s.get('https://string-db.org/cgi/download.pl')
        absolute_links = r.html.absolute_links
        string_dir = os.path.join(self.installation_directory, 'data')

        if self.string_links == 'download':
            self.tell('Downloading STRING interactions database')
            string_url = None
            for l in absolute_links:
                if 'protein.links.full' in l:
                    string_url = l
                    break
            if string_url:
                command = "wget '{url}' -P '{dir}'".format(url=string_url, dir=string_dir)
                subprocess.call(command, shell=True)
                self.string_links = os.path.join(string_dir, os.path.basename(string_url))
        else:
            self.tell('STRING interactions file provided, adding', self.string_links, ' to configuration file')
        Configuration.CONFIG.set('databases', 'string_links', self.string_links)

        if self.string_sequences == 'download':
            self.tell('Downloading STRING sequences database')
            string_url = None
            for l in absolute_links:
                if 'protein.sequences' in l:
                    string_url = l
                    break
            command = "wget '{url}' -P '{dir}'".format(url=string_url, dir=string_dir)
            subprocess.call(command, shell=True)
            self.string_sequences = os.path.join(string_dir, os.path.basename(string_url))
        else:
            self.tell('STRING sequences file provided, adding', self.string_sequences, ' to configuration file')
        Configuration.CONFIG.set('databases', 'string_sequences', self.string_sequences)

        if self.string_species == 'download':
            self.tell('Downloading STRING sequences database')
            string_url = None
            for l in absolute_links:
                if 'species.v' in l:
                    string_url = l
                    break
            command = "wget '{url}' -P '{dir}'".format(url=string_url, dir=string_dir)
            subprocess.call(command, shell=True)
            self.string_species = os.path.join(string_dir, os.path.basename(string_url))
        else:
            self.tell('STRING species file provided, adding', self.string_species, ' to configuration file')
        Configuration.CONFIG.set('databases', 'string_species', self.string_species)

        uniprot_dir = os.path.join(self.installation_directory, 'data/UniprotKB')
        if self.uniprot_sprot == 'download':
            self.tell('Downloading UniProt SwissProt file')
            uniprot_url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/' \
                          'uniprot_sprot.fasta.gz'
            command = "wget '{url}' -P '{dir}'".format(url=uniprot_url, dir=uniprot_dir)
            subprocess.call(command, shell=True)
            self.uniprot_sprot = os.path.join(uniprot_dir, 'uniprot_sprot.fasta.gz')
        else:
            self.tell('UniProt SwissProt file provided, adding', self.uniprot_sprot, ' to configuration file')
        Configuration.CONFIG.set('databases', 'uniprot_sprot', self.uniprot_sprot)

        if self.uniprot_goa == 'download':
            self.tell('Downloading UniProt GOA file')
            uniprot_url = 'ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz'
            command = "wget '{url}' -P '{dir}'".format(url=uniprot_url, dir=uniprot_dir)
            subprocess.call(command, shell=True)
            self.uniprot_goa = os.path.join(uniprot_dir, 'goa_uniprot_all.gaf.gz')
        else:
            self.tell('UniProt GOA file provided, adding', self.uniprot_goa, ' to configuration file')
        Configuration.CONFIG.set('databases', 'uniprot_goa', self.uniprot_goa)

    def process_files(self):
        # look for the uncompressed version of every of the downloaded datasets
        self.tell('Uncompressing STRING')
        # STRING interactions
        uncompressed_file = self.string_links.split('.gz')[0]
        if not os.path.exists(uncompressed_file):
            command = 'gunzip ' + self.string_links
            subprocess.call(command, shell=True)
        else:
            self.tell('Found uncompressed STRING interactions file, skipping')

        core_ids_file = os.path.join(self.installation_directory, 'data/coreIds')
        core_ids = []
        if os.path.exists(core_ids_file):
            self.tell('Found coreIds file, skipping inference')
        else:
            self.tell('Extracting core IDs')
            for line in open(self.string_species, 'r'):
                if line[0] == '#':
                    continue
                fields = line.split('\t')
                if fields[1] == 'core':
                    core_ids.append(fields[0])
            f = open(core_ids_file, 'w')
            f.write('\n'.join(core_ids))
            f.close()
        self.tell('Found', len(core_ids), 'core IDs')

        # STRING sequences
        uncompressed_file = self.string_sequences.split('.gz')[0]
        if not os.path.exists(uncompressed_file):
            command = 'gunzip ' + self.string_sequences
            subprocess.call(command, shell=True)
        else:
            self.tell('Found uncompressed STRING sequences file, skipping')

        self.tell('Partitioning STRING sequences file')
        current_organism = ''
        fp = None
        with open(uncompressed_file, 'r') as string_sequences:
            for line in string_sequences:
                if line.startswith('>'):
                    taxid = line[1:].split('.')[0]
                    if current_organism != taxid:  # different organism
                        if fp:
                            fp.close()
                            fp = None
                        if taxid in core_ids:
                            current_organism = taxid
                            fn = os.path.join(self.installation_directory, 'data/STRINGSequences/'+taxid+'.faa')
                            if not os.path.exists(fn):
                                fp = open(fn, 'a')
                                fp.write(line)
                    else:
                        if fp:
                            fp.write(line)
                else:
                    if fp:
                        fp.write(line)
        if fp:
            fp.close()

        self.tell('Creating BLAST databases from STRING sequences')
        seq_dir = os.path.join(self.installation_directory, 'data/STRINGSequences')
        for fasta in os.listdir(seq_dir):
            if re.match('.*\.faa$', fasta):
                f = os.path.abspath(os.path.join(seq_dir, fasta))
                if not os.path.exists(f+'.phr'):
                    subprocess.call('makeblastdb -in {fasta} -out {fasta} -dbtype prot'.format(fasta=f), shell=True)

        self.tell('Uncompressing UniProtKB')
        # UniProt SwissProt
        uncompressed_file = self.uniprot_sprot.split('.gz')[0]
        if not os.path.exists(uncompressed_file):
            command = 'gunzip ' + self.uniprot_sprot
            subprocess.call(command, shell=True)
        else:
            self.tell('Found uncompressed UniProt SwissProt file, skipping')

        # UniProt GOA
        uncompressed_file = self.uniprot_goa.split('.gz')[0]
        if not os.path.exists(uncompressed_file):
            command = 'gunzip ' + self.uniprot_goa
            subprocess.call(command, shell=True)
        else:
            self.tell('Found uncompressed UniProt GOA file, skipping')

        # Filter UniProt GOA
        if self.filtered_goa == 'infer':
            self.tell('filtering UniProt')
            self.filter_uniprot()

    def filter_uniprot(self):
        self.tell('Gathering UniProt SwissProt IDs')
        ids = []
        proteins_with_function = set()
        for line in open(self.uniprot_sprot.split('.gz')[0], 'r'):
            if line.startswith('>'):
                ids.append(line.split('|')[1])
        self.tell('Gathered', len(ids), 'SwissProt IDs')

        self.tell('Filtering UniProt GOA')
        experimental_goa = self.uniprot_goa.split('.gz')[0] + '.exp'

        num_lines = Utilities.wccount(experimental_goa)

        if self.evidence_codes == 'experimental':
            self.evidence_codes = GeneOntology.GeneOntology.EXPERIMENTAL_EVIDENCE_CODES
        Configuration.CONFIG.set('options', 'evidence_codes', ','.join(self.evidence_codes))

        self.tell('Accepting the following evindece codes:', self.evidence_codes)

        # filter evidence codes using awk
        command = "awk '$6~/" + '|'.join(self.evidence_codes)+"/{print $0}' " + \
                  self.uniprot_goa.split('.gz')[0] + " > " + experimental_goa
        subprocess.call(command, shell=True)
        # keep only SwissProt annotations
        self.filtered_goa = os.path.join(self.installation_directory, 'data/UniprotKB/filtered_goa')
        fg = open(self.filtered_goa, 'w')
        i = 0
        for line in open(experimental_goa, 'r'):
            i+=1
            if i % 10000 == 0:
                self.tell('Process:', i/num_lines*100.0,'%')
            if not line.startswith('!'):
                fields = line.split('\t')
                if fields[1] in ids:
                    proteins_with_function.add(fields[1])
                    fg.write(line)
        fg.close()
        Configuration.CONFIG.set('databases', 'filtered_goa', self.filtered_goa)

        self.tell('Filtering UniProt SwissProt')
        self.filtered_sprot = os.path.join(self.installation_directory, 'data/UniprotKB/filtered_sprot')
        fs = open(self.filtered_sprot, 'w')
        writing = False
        for line in open(self.uniprot_sprot.split('.gz')[0], 'r'):
            if line.startswith('>'):
                prot_id = line.split('|')[1]
                writing = prot_id in proteins_with_function
                if writing:
                    fs.write(line)
            elif writing:
                fs.write(line)
        fs.close()
        Configuration.CONFIG.set('databases', 'filtered_sprot', self.filtered_sprot)

    def summary_and_continue(self):
        summary = 'These are the loaded values \n'
        summary += '\tInstallation directory:\t\t' + self.installation_directory + '\n\n'

        summary += '\tPath to InterPro executable:\t' + self.interpro + '\n'
        summary += '\tPath to HMMer:\t\t\t' + self.hmmer + '\n'
        summary += '\tPath to blastp:\t\t\t' + self.blastp + '\n'
        summary += '\tPath to makeblastdb:\t\t' + self.makeblastdb + '\n\n'

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

        summary += '\tEvidence Codes:\t\t\t' + str(self.evidence_codes) + '\n'

        summary += 'Do you want to continue with these settings?'

        if Utilities.query_yes_no(summary):
            return
        sys.exit()
