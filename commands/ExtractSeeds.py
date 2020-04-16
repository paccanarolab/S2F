import pandas as pd
from Utils import FancyApp
from diffusion import Diffusion
from scipy import sparse
import os


class ExtractSeeds(FancyApp.FancyApp):

    def __init__(self, args):
        super(ExtractSeeds, self).__init__()
        self.prediction_directory = args.prediction_directory
        self.alias = self.alias

    def run(self):
        self.tell('Reading proteins index')
        proteins = pd.read_pickle(os.path.join(self.prediction_directory,
                                               'proteins.df'))
        self.tell('Reading terms index')
        terms = pd.read_pickle(os.path.join(self.prediction_directory,
                                            'terms.df'))
        self.tell('Loading interpro seed file')
        interpro = sparse.load_npz(os.path.join(self.prediction_directory,
                                                '../../seeds/interpro',
                                                f'{self.alias}.npz'))
        self.tell('Loading hmmer seed file')
        hmmer = sparse.load_npz(os.path.join(self.prediction_directory,
                                             '../../seeds/hmmer',
                                             f'{self.alias}.npz'))
        self.tell('saving text version of interpro seed')
        Diffusion._write_results(interpro, proteins, terms,
                                 os.path.join(self.prediction_directory,
                                              '../../seeds/interpro',
                                              f'{self.alias}.seed.txt'))
        self.tell('saving text version of hmmer seed')
        Diffusion._write_results(hmmer, proteins, terms,
                                 os.path.join(self.prediction_directory,
                                              '../../seeds/interpro',
                                              f'{self.alias}.seed.txt'))
        self.tell('done')
