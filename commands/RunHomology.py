import os
import subprocess

from graphs import homology
from Utils import ColourClass, FancyApp, Utilities


class RunHomology(FancyApp.FancyApp):

    def __init__(self, args):
        super(RunHomology, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_GREEN
        self.fasta = args.fasta
        self.output_dir = args.output_dir
        self.alias = args.alias
        self.tell('extracting proteins from fasta')
        self.proteins = Utilities.extract_indices_from_fasta(self.fasta)

    def run(self):
        subprocess.call(f"makeblastdb -in {self.fasta} -out {self.fasta} " +
                        "-dbtype prot",
                        shell=True)
        h = homology.Homology(
            self.fasta, self.proteins, self.output_dir, self.alias)
        self.tell('computing homology graph')
        h.compute_graph()
        self.tell('writing homology graph file')
        h.write_graph(os.path.join(self.output_dir, f'{self.alias}.homology'))
