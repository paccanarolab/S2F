from GOTool import GeneOntology
from Utils import FancyApp, ProgressBar
import pandas as pd
import numpy as np
import os

class RescoreContinuous(FancyApp.FancyApp):

    def __init__(self, args):
        super(RescoreContinuous, self).__init__()
        self.obo = args.obo
        self.prediction = os.path.expanduser(args.prediction)
        self.outdir = os.path.expanduser(args.outdir)
        self.th_start = args.th_start
        self.go = None

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
        self.go = GeneOntology.GeneOntology(self.obo, verbose=True)
        self.go.build_structure()
        self.tell('Loading prediction file...')
        prediction_df = pd.read_csv(self.prediction, header=None, sep='\t',
                                    names=['protein', 'go_term', 'score'],
                                    dtype={'protein':str, 'go_term':str, 'score':np.float32})
        filtered_df = None
        proteins = prediction_df.protein.unique()
        if self.__verbose__:
            prog = ProgressBar.ProgressBar(0, len(proteins), 77, mode='dynamic', char='-')
        for protein in proteins:
            predictions = prediction_df[prediction_df.protein == protein]
            th = self.th_start
            while len(predictions) > 1500:
                predictions = self.filter_predictions(predictions, th, protein)
            if filtered_df is None:
                filtered_df = predictions
            else:
                filtered_df = filtered_df.append(predictions, ignore_index=True)
            if self.__verbose__:
                prog.increment_amount()
                prog.print_bar()
        if self.__verbose__:
            prog.finish_bar()
        
        outfile = os.join(self.outdir, os.path.basename(self.prediction) + '_filtered')
        self.tell('Saving prediction to file...')
        filtered_df.to_csv(outfile, sep='\t', header=None)




