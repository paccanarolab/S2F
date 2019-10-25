from Utils import *

from GOTool import GeneOntology

import pandas as pd
import os


class BuildClamp(FancyApp.FancyApp):

    def __init__(self, args):
        super(BuildClamp, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_GREEN
        # this file requires no configuration, only paths
        self.config_file = os.path.expanduser(args.config_file)
        Configuration.load_configuration(self.config_file)
        self.fasta = os.path.expanduser(args.fasta)
        if args.output == 'infer':
            self.output = self.fasta + '.clamp'
        else:
            self.output = os.path.expanduser(args.output)
        self.proteins = Utilities.extract_indices_from_fasta(self.fasta, Utilities.keep_uniprot_accession)
        self.evidence_codes = args.evidence_codes

    def run(self):
        uniprot_goa = Configuration.CONFIG.get('databases', 'uniprot_goa')
        uniprot_goa = uniprot_goa.split('.gz')[0]
        if self.evidence_codes == 'experimental':  # in case the user request experimental, we cached that
            uniprot_goa = uniprot_goa + '.exp'
            self.evidence_codes = GeneOntology.GeneOntology.EXPERIMENTAL_EVIDENCE_CODES
        self.tell('Filtering GOA file')
        out = open(self.output, 'w')
        num_lines = Utilities.line_count(uniprot_goa)
        if self.__verbose__:
            prog = ProgressBar.ProgressBar(0, num_lines, 77, mode='dynamic', char='-')
        for line in open(uniprot_goa):
            fields = line.split('\t')
            if fields[1] in self.proteins.index and fields[6] in self.evidence_codes:
                out.write(line)
            if self.__verbose__:
                prog.increment_amount()
                prog.print_bar()
        if self.__verbose__:
            prog.finish_bar()
        out.close()
        self.tell('done, written output file:', self.output)


