import pandas as pd
from GOTool import GeneOntology
from seeds import interpro
from Utils import FancyApp


class SeedFromInterPro(FancyApp.FancyApp):

    def __init__(self, args):
        super(SeedFromInterPro, self).__init__()
        self.obo = args.obo
        self.interpro_file = args.interpro_file
        self.output = args.output

    def run(self):
        self.tell('Building Gene Ontology')
        go = GeneOntology.GeneOntology(self.obo, verbose=True)
        go.build_structure()
        terms = pd.DataFrame(list(enumerate(sorted(go.terms.keys()))))
        self.tell('Extracting protein set from InterPro file')
        proteins = set()
        with open(self.interpro_file) as ip_file:
            for line in ip_file:
                proteins.add(line.split('\t')[0])
        proteins = pd.DataFrame(list(enumerate(sorted(list(proteins)))))
        self.tell('Parsing the seed file...')
        seeder = interpro.InterProSeed(self.interpro_file,
                                       proteins,
                                       terms,
                                       go)
        seeder.process_output()
        assignment = seeder.get_seed(return_pandas_assignment=True)
        self.tell('Saving seed file into', self.output, '...')
        assignment.to_csv(self.output, sep='\t', header=None, index=False)
        self.tell('done')
