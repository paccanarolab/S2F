from Utils import ColourClass, FancyApp, ProgressBar
from GOTool import AnnotationParser

import datetime
import os

__author__ = 'Mateo Torres'
__email__ = 'Mateo.Torres.2015@live.rhul.ac.uk'
__copyright__ = 'Copyright (c) 2020, Mateo Torres'
__license__ = 'MIT'
__version__ = '0.1'


class GOAHistory(FancyApp):

    def __init__(self, args):
        super(GOAHistory, self).__init__()
        self.colour = ColourClass.bcolors.MAGENTA
        self.now = datetime.datetime.now()
        self.year = int(args.get('year', self.now.year))
        self.annotation_file = args.annotation_file
        self.output_dir = args.output_dir
        self.organism = args.organism

    def extract_annotations(self):
        parser = AnnotationParser.AnnotationFile(self.annotation_file,
                                                 self.organism)
        self.tell('Filtering annotations, keeping only annotations before',
                  self.year)
        counter = 0
        fname = os.path.join(self.output_dir, self.organism + '.gaf.' +
                             str(self.year))
        outfile = open(fname, 'w')
        for annotation in parser:
            date = datetime.datetime.strptime(annotation.date, '%Y%m%d')
            year = date.year
            if year >= self.year:
                if len(annotation.qualifiers) > 0 and \
                   annotation.qualifiers[0] != 'NOT':
                    continue
                else:
                    outfile.write(f'{annotation.go_id}\t'
                                  f'{annotation.db_object_id}\n')
        outfile.close()



