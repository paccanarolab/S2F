from GOTool import GeneOntology
from Utils import FancyApp
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

    def distribute_scores(self, df):
        pass

    def filter_predictions(self, df, th):
        res = pd.DataFrame()
        res.to_csv('asdf.txt',)
        return res

    def run(self):
        self.tell('Building Gene Ontology')
        go = GeneOntology.GeneOntology(self.obo, verbose=True)
        go.build_structure()
        self.tell('Loading prediction file...')
        prediction_df = pd.read_csv(self.prediction, header=None, sep='\t'
                                    names=['protein', 'go_term', 'score'],
                                    dtype={'protein':str, 'go_term':str, 'score':np.float32})
        filtered_df = None
        proteins = prediction_df.proteins.unique()
        for protein in proteins:
            predictions = prediction_df[prediction_df.protein == protein]
            th = self.th_start
            while len(predictions) > 1500:
                predictions = self.filter_predictions(predictions, th)
            predictions = self.distribute_scores(predictions)
            if filtered_df is None:
                filtered_df = predictions
            else:
                filtered_df = filtered_df.append(predictions, ignore_index=True)

        outfile = os.join(self.outdir, os.path.basename(self.prediction) + '_filtered')
        self.tell('Saving prediction to file...')
        filtered_df.to_csv(outfile, sep='\t', header=None)




