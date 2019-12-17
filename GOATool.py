"""
GOATool - Main Script
"""
from GOTool.GOAHistory import GOAHistory
import argparse
import os

aparser = argparse.ArgumentParser()
required_arguments = aparser.add_argument_group('required arguments')
required_arguments.add_argument('--output-dir', '--od',
                                help='path where the filtered gaf will be '
                                     'written',required=True)
required_arguments.add_argument('--annotation_file', '--gaf',
                                help='GOA annotation file in GAF format',
                                required=True)
required_arguments.add_argument('--organism',
                                help='Organism name', required=True)
required_arguments.add_argument('--year', required=True)


args = aparser.parse_args()

hg = GOAHistory(args)
hg.extract_annotations()