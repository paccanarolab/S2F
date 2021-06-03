from Measures import BuildMatrices, measures
from Utils import FancyApp, ColourClass
from GOTool import GeneOntology

from collections import defaultdict
from scipy import sparse
import numpy as np
import pandas as pd

import argparse
import os
import matplotlib
matplotlib.use('Agg')


def save_list_to_txt(l, file_name):
    with open(file_name, "w") as out:
        out.write("\n".join(map(str, l)))


def tell(*args, **kwargs):
    FancyApp.FancyApp.yell(
        ColourClass.bcolors.BOLD_GREEN, 'evaluator', *args, **kwargs)


aparser = argparse.ArgumentParser()
required_arguments = aparser.add_argument_group('required arguments')
required_arguments.add_argument('--prediction',
                                help='prediction file', required=True)
required_arguments.add_argument('--goa', help='goa file',
                                required=True)
required_arguments.add_argument('--proteins-idx',
                                help='path to the proteins index file, used'
                                'only when prediction format is numpy',
                                required=False)
required_arguments.add_argument('--terms-idx',
                                help='path to the terms index file, used'
                                'only when prediction format is npz',
                                required=False)
required_arguments.add_argument('--obo',
                                help='go.obo file', required=True)
required_arguments.add_argument('--outputdir',
                                help='output directory', required=True)
required_arguments.add_argument('--output-matrices', action='store_true')
required_arguments.add_argument('--compute-metrics', action='store_true')
required_arguments.add_argument('--prediction-fmt',
                                help='how the prediction file should be read',
                                default='txt',
                                choices=['txt', 'npz'])


args = aparser.parse_args()


# we have these values per threshold, so they are unsafe for plotting
unsafe_metrics = ['roc', 'pr', 'tp', 'fp', 'tn', 'fn', 's', 'ru', 'mi',
                  'AUC per-gene raw', 'AUPR per-gene raw',
                  'Precision at 0.2 Recall per-gene raw', 'F_max per-gene raw',
                  'NDCG per-gene raw', 'Jaccard per-gene raw',
                  'smin per-gene raw', 'AUC per-term raw', 'AUPR per-term raw',
                  'Precision at 0.2 Recall per-term raw', 'F_max per-term raw',
                  'NDCG per-term raw', 'Jaccard per-term raw',
                  'smin per-term raw']

tell('Loading Gene Ontology and annotations')
go = GeneOntology.GeneOntology(args.obo, verbose=True)
go.build_structure()
go.load_annotation_file(args.goa, 'GOA')
go.up_propagate_annotations('GOA')
annotations = go.get_annotations('GOA')

if args.prediction_fmt == 'txt':
    proteins = set()
    terms = set()
    tell('Inferring proteins and terms from prediction file')
    for line in open(args.prediction):
        fields = line.strip().split()
        proteins.add(fields[0])
        terms.add(fields[1])

    tell('Inferring proteins and terms from GOA')
    proteins |= set(annotations['Protein'].unique())
    terms |= set(annotations['GO ID'].unique())

    tell('Sorting terms...')
    terms = sorted(list(terms))
    tell('Indexing terms...')
    terms_idx = {term: i for i, term in enumerate(terms)}

    tell('Sorting proteins...')
    proteins = sorted(list(proteins))
    tell('Indexing proteins...')
    proteins_idx = {protein: i for i, protein in enumerate(proteins)}
else:
    tell('Loading protein index')
    proteins_idx = pd.read_pickle(args.proteins_idx).to_dict()['protein idx']
    proteins = set(proteins_idx.keys())
    tell('loading terms index')
    terms_idx = pd.read_pickle(args.terms_idx).to_dict()['term idx']
    terms = set(terms_idx.keys())

prediction = np.zeros((len(proteins), len(terms)))
goa = np.zeros((len(proteins), len(terms)))
information_content = np.zeros(len(terms))
if args.prediction_fmt == 'text':
    tell('Building prediction matrix...')
    for line in open(args.prediction):
        fields = line.strip().split()
        score = float(fields[2])
        if score > 0:
            prediction[proteins_idx[fields[0]],
                       terms_idx[fields[1]]] = score
else:
    tell('Loading prediction matrix...')
    prediction = sparse.load_npz(args.prediction).toarray()
tell('Building GOA matrix')
for i, row in annotations.iterrows():
    goa[proteins_idx[row['Protein']],
        terms_idx[row['GO ID']]] = 1.0

if not os.path.exists(args.outputdir):
    os.makedirs(args.outputdir)

if args.output_matrices:
    tell('saving prediction matrix')
    np.save(os.path.join(args.outputdir,'prediction.npy'), prediction)

    tell('saving GOA matrix')
    np.save(os.path.join(args.outputdir,'goa.npy'), goa)

if args.compute_metrics:
    tell('Calculating Information Content')
    for term, i in terms_idx.items():
        information_content[i] = go.find_term(term).information_content('GOA')

    tell('saving indices')
    prots_df = pd.DataFrame({'protein':list(proteins_idx.keys()), 
                            'index':list(proteins_idx.values())})
    prots_df.to_csv(os.path.join(args.outputdir, 'proteins.tsv'), sep='\t')
    terms_df = pd.DataFrame({'term':list(terms_idx.keys()),
                            'index':list(terms_idx.values())})
    terms_df.to_csv(os.path.join(args.outputdir, 'terms.tsv'), sep='\t')

    low, high = 1, 1000000

    sumcol = np.sum(goa, axis=0)
    sumrow = np.sum(goa, axis=1)
    tell('Extracting gold standard')

    pred = prediction[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)]
    if len(np.unique(pred)) > 10000:
        pred = np.around(pred, decimals=4)
    gold_standard = goa[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)]
    ic = information_content[(sumcol >= low) & (sumcol <= high)]

    tell('gold_standard', gold_standard.shape)
    tell('information_content', ic.shape)

    data = defaultdict(list)
    meas = measures.HX_py(pred, ic, 1234)
    tell('computing OVERALL')
    res = meas.compute_overall(gold_standard)
    for m in res.keys():
        data[m].append(res[m])

    tell('computing per-gene')
    res = meas.compute_per_gene(gold_standard)
    for m in res.keys():
        data[m].append(res[m])

    tell('computing per-term')
    res = meas.compute_per_term(gold_standard)
    for m in res.keys():
        data[m].append(res[m])

    tell('saving metrics')
    # save all metrics
    metrics = []
    for metric, values in data.items():
        metrics.append((metric, values))
    metrics_df = pd.DataFrame(
        metrics, columns=['metric', 'value'])
    metrics_df.to_pickle(os.path.join(args.outputdir,'metrics_df.pkl'))
