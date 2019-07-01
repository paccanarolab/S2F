from Utils import *
from GOTool import GeneOntology

from graphs import Graph
from graphs import combination

from scipy import sparse
from scipy.sparse.linalg import lsqr
from sklearn.preprocessing import binarize

import pandas as pd
import numpy as np
import os


class CombineSeeds(FancyApp.FancyApp):

    def __init__(self, args):
        super(CombineSeeds, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_GREEN
        self.seed_files = args.seed_files
        self.seed_separator = args.seed_separator
        self.output = args.output
        self.coefficients = args.coefficients
        if self.coefficients == 'infer':
            n = len(self.seed_files)
            self.coefficients = [1/n] * n
        elif len(self.coefficients) != len(self.seed_files):
            raise ValueError("The number of coefficients is different than the number of seeds")

    def run(self):
        seeds = []
        names = ['protein', 'goterm', 'score']
        self.tell('Reading seed files')
        for seed_file in self.seed_files:
            seeds.append(pd.read_csv(seed_file, sep=self.seed_separator, names=names))
        seeds[0]['score'] = self.coefficients[0] * seeds[0]['score']
        aggregated_seed = seeds[0]
        for i, seed in enumerate(seeds[1:]):
            aggregated_seed = aggregated_seed.merge(seed,
                                                    on=['protein', 'goterm'],
                                                    how='outer',
                                                    suffixes=['_1', '_2']).fillna(0)
            aggregated_seed['score'] = aggregated_seed['score_1'] + self.coefficients[i+1] * aggregated_seed['score_2']
            aggregated_seed.drop(columns=['score_1', 'score_2'])
        self.tell('writing aggregated seed')
        aggregated_seed.to_csv(self.output, sep=self.seed_separator, header=False)