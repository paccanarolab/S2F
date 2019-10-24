from GOTool import GeneOntology
from seeds import hmmer
from Utils import FancyApp


class SeedFromHMMER(FancyApp.FancyApp):

    def __init__(self, args):
        super(SeedFromHMMER, self).__init__()
        self.obo = args.obo
        self.evalue_file = args.evalue_file
        self.output = args.output
        self.threshold = args.threshold

    def run(self):
        self.tell('Building Gene Ontology')
        go = GeneOntology.GeneOntology(self.obo, verbose=True)
        go.build_structure()
        self.tell('Parsing the seed file...')
        seeder = hmmer.HMMerSeed('', None, None, go, [], '')
        seeder.process_output(evalue_file=self.evalue_file)
        assignment = seeder.get_seed(seed_threshhold=self.threshold,
                                     return_pandas_assignment=True)
        self.tell('Saving seed file into', self.output, '...')
        assignment.to_csv(self.output, sep='\t', header=None, index=False)
        self.tell('done')
