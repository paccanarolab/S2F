import pandas as pd
from Utils import ColourClass, FancyApp

class TopPredictions(FancyApp.FancyApp):
    
    def __init__(self, args):
        super(TopPredictions, self).__init__()
        self.colour = ColourClass.bcolors.BOLD_GREEN
        self.prediction_file = args.prediction_file
        self.n = args.N
        self.output = args.output
    
    def run(self):
        self.tell(f'Loading prediction file {self.prediction_file}')
        prediction = pd.read_csv(self.prediction_file,
                                 sep='\t',
                                 names=['protein', 'goterm', 'score'])
        self.tell(f'getting top {self.n} predictions per-protein')
        top_n = (prediction
                 .sort_values('score', ascending=False)
                 .groupby('protein').head(self.n).reset_index())
        self.tell(f'Saving results into {self.output}')
        top_n.to_csv(self.output, sep='\t', index=False, header=False)
