from GOTool import GeneOntology
from Utils import FancyApp, ProgressBar
from multiprocessing import Pool

import pandas as pd
import numpy as np
import os

PREDICTION_DF = None
GO = None

def filter_predictions_worker(th, protein):
    df = PREDICTION_DF[PREDICTION_DF.protein == protein]
    data = {'protein':[], 'go_term':[], 'score':[]}
    term_cache = set()
    for go_term in df.go_term.unique():
        term = GO.find_term(go_term)
        score = df[df.go_term == go_term].at[0, 'score']
        children = term.get_children()
        keep = True
        for c in children:
            if len(df[df.go_term == c.name]) < 1:
                continue
            c_score = df[df.go_term == c.name].at[0, 'score']
            if c_score <= score + th:
                keep = False
                break
        if keep:
            data['protein'].append(protein)
            data['go_term'].append(go_term)
            data['score'].append(score)
    res = pd.DataFrame(data)
    return res

class RescoreContinuous(FancyApp.FancyApp):

    def __init__(self, args):
        super(RescoreContinuous, self).__init__()
        self.obo = args.obo
        self.prediction = os.path.expanduser(args.prediction)
        self.outdir = os.path.expanduser(args.outdir)
        self.th_start = args.th_start
        self.go = None
        self.cpu = len(os.sched_getaffinity(0))

    def filter_predictions(self, df, th, protein):
        df.to_csv('test.df', sep='\t', header=None)
        data = {'protein':[], 'go_term':[], 'score':[]}
        term_cache = set()
        for go_term in df.go_term.unique():
            term = self.go.find_term(go_term)
            score = df[df.go_term == go_term].at[0, 'score']
            children = term.get_children()
            keep = True
            for c in children:
                if len(df[df.go_term == c.name]) < 1:
                    continue
                c_score = df[df.go_term == c.name].at[0, 'score']
                if c_score <= score + th:
                    keep = False
                    break
            if keep:
                data['protein'].append(protein)
                data['go_term'].append(go_term)
                data['score'].append(score)
        res = pd.DataFrame(data)
        return res

    def run(self):
        self.tell('Building Gene Ontology')
        GO = GeneOntology.GeneOntology(self.obo, verbose=True)
        GO.build_structure()
        self.tell('Loading prediction file...')
        PREDICTION_DF = pd.read_csv(self.prediction, header=None, sep='\t',
                                    names=['protein', 'go_term', 'score'],
                                    dtype={'protein':str, 'go_term':str, 'score':np.float32})
        filtered_df = None
        proteins = PREDICTION_DF.protein.unique()
        params = []
        self.tell('gathering parameters for parallel pool')
        if self.__verbose__:
            prog = ProgressBar.ProgressBar(0, len(proteins), 77, mode='dynamic', char='-')
        for protein in proteins:
            th = self.th_start
            params.append([th, protein])
            if self.__verbose__:
                prog.increment_amount()
                prog.print_bar()
        if self.__verbose__:
            prog.finish_bar()
        self.tell(f'Filtering DataFrames using {self.cpu} cores')
        with Pool(self.cpu) as p:
            p.starmap(params, filter_predictions_worker)
        
        outfile = os.join(self.outdir, os.path.basename(self.prediction) + '_filtered')
        self.tell('Saving prediction to file...')
        filtered_df.to_csv(outfile, sep='\t', header=None)




