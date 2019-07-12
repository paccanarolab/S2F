from Measures import BuildMatrices, measures
from Utils import FancyApp, ColourClass

from collections import defaultdict

import numpy as np
import pandas as pd

import argparse
import os
import re
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
# from Utils.notification import load_configuration, get_api


def tell(*args, **kwargs):
    FancyApp.FancyApp.yell(ColourClass.bcolors.BOLD_GREEN, 'plotter', *args, **kwargs)


aparser = argparse.ArgumentParser()
required_arguments = aparser.add_argument_group('required arguments')
required_arguments.add_argument('--outputs-dir', '--od',
                                help='tsv that links tax ids to output folders, not the full path'
                                     'but only the name of the containig directory to avoid confussion', required=True)
required_arguments.add_argument('--cm-diffusion-dir', '--cm',
                                help='directory containing the results computed by the '
                                     'consistency method diffusion', required=True)
required_arguments.add_argument('--plot-individual', action='store_true')
required_arguments.add_argument('--compute-metrics', action='store_true')


args = aparser.parse_args()

output_dirs = pd.read_csv(args.outputs_dir, sep='\t',
                          names=['Tax ID', 'output'],
                          dtype={'Tax ID': 'str', 'output': 'str'})

# This is a hack, and should really be an argument. This is a plot generation script though, so no much thought into it.
if os.name == 'nt':
    dataroot = r'A:\COMMON\PROJECTS\S2F\data'
else:
    dataroot = '/home/paccanaro/COMMON/PROJECTS/S2F/data/'

if os.name == 'nt':
    DATA_DIRECTORY = r'G:\PaccanaroLab\s2f\data\bacteria_selection/'
else:
    DATA_DIRECTORY = '../../../data/bacteria_selection/'

obo = os.path.join(dataroot, 'INPUT/go.obo')
OUTPUT_BASE = os.path.join(dataroot, 'output/')
GOA_DIRECTORY = os.path.join(DATA_DIRECTORY, 'selected_goa/')
FASTA_DIRECTORY = os.path.join(DATA_DIRECTORY, 'selected_fasta/')
PLOTS_DIRECTORY = os.path.join(DATA_DIRECTORY, 'comparison_plots/')
MATRICES_DIRECTORY = os.path.join(PLOTS_DIRECTORY, 'matrices/')

# load the selected organisms pandas fil
final_selection = pd.read_pickle(DATA_DIRECTORY+'final_selection.pkl')

all_organisms = pd.DataFrame()

# for coloring the plots
flatui = ["#705319", "#e5aa34", "#911f76", "#e534bb", "#1c7c69", "#31d6b5", ]
flatui3 = ["#e5aa34", "#e534bb", "#31d6b5", ]

# for nice labelling
labels = [
    ('I+H+D', 'S2F'),
    ('I+H+CM', 'InterPro & HMMER + CM'),
    ('I+H', 'InterPro + HMMER'),
    ('I', 'InterPro'),
    ('I+D', 'InterPro w/ diffusion'),
    ('I+CM', 'InterPro + CM'),
    ('H', 'HMMER'),
    ('H+D', 'HMMER w/ diffusion'),
    ('H+CM', 'HMMER + CM'),
]

# we have these values per threshold, so they are unsafe for plotting
unsafe_metrics = ['s', 'ru', 'mi', 'roc', 'pr']

condition = final_selection['Tax ID'].isin(output_dirs['Tax ID'])

for i, organism in final_selection[condition].sort_values(by='Tax ID').iterrows():
    tell('processing organism: ', organism['Organism'], '(', organism['Tax ID'], ')...')
    if organism['Tax ID'] in ['1392', '243277', '272942']:  # we exclude this organism because it has too many proteins for us to run it
        tell('skipping due to unfinished computation')
        continue

    interpro_file = os.path.join(dataroot,'seeds/interpro/S2F-2018-'+organism['Tax ID']+'.matrix')
    hmmer_file = os.path.join(dataroot, 'seeds/hmmer/S2F-2018-'+organism['Tax ID']+'.matrix')
    # the output directory will be the one with the latest date in the directory name
    output_dir = os.path.join(OUTPUT_BASE, output_dirs[output_dirs['Tax ID'] == organism['Tax ID']]['output'].values[0])
    tell('loading output_files from', output_dir)
    interpro_diff_file = os.path.join(output_dir, 'interpro.diffusion')
    hmmer_diff_file = os.path.join(output_dir, 'hmmer.diffusion')
    hmmer_cm_file = os.path.join(args.cm_diffusion_dir, 'diffusion-CM-hmmer-' + str(organism['Tax ID']))
    interpro_cm_file = os.path.join(args.cm_diffusion_dir, 'diffusion-CM-' + str(organism['Tax ID']))
    s2f_file = os.path.join(output_dir, 'final.diffusion')
    goa = os.path.join(GOA_DIRECTORY, organism['File'])
    fasta = os.path.join(FASTA_DIRECTORY, organism['Tax ID'] + '.fasta')

    # load or build matrices
    org_dir = os.path.join(MATRICES_DIRECTORY, organism['Tax ID'] + '/')

    if os.path.exists(org_dir + 'metrics_df.pkl'):
        tell('Found saved metrics file, loading...')
        metrics_df = pd.read_pickle(org_dir + 'metrics_df.pkl')
    else:
        if os.path.exists(org_dir):
            tell('Found precomputed matrices, loading...')
            goa_values = np.load(org_dir + 'goa_values.npy')
            hmmer_diff_values = np.load(org_dir + 'hmmer_diff_values.npy')
            hmmer_cm_values = np.load(org_dir + 'hmmer_cm_values.npy')
            hmmer_values = np.load(org_dir + 'hmmer_values.npy')
            interpro_diff_values = np.load(org_dir + 'interpro_diff_values.npy')
            interpro_cm_values = np.load(org_dir + 'interpro_cm_values.npy')
            interpro_values = np.load(org_dir + 'interpro_values.npy')
            s2f_values = np.load(org_dir + 's2f_values.npy')
            information_content = np.load(org_dir + 'information_content.npy')
        else:
            os.makedirs(org_dir)
            matrix_builder = BuildMatrices.BuildMatrices(
                interpro_file, interpro_diff_file, interpro_cm_file,
                hmmer_file, hmmer_diff_file, hmmer_cm_file, s2f_file, goa, fasta, obo)
            np.save(org_dir + 'goa_values.npy', matrix_builder.goa_values)
            np.save(org_dir + 'hmmer_diff_values.npy', matrix_builder.hmmer_diff_values)
            np.save(org_dir + 'hmmer_cm_values.npy', matrix_builder.hmmer_cm_values)
            np.save(org_dir + 'hmmer_values.npy', matrix_builder.hmmer_values)
            np.save(org_dir + 'interpro_diff_values.npy', matrix_builder.interpro_diff_values)
            np.save(org_dir + 'interpro_cm_values.npy', matrix_builder.interpro_cm_values)
            np.save(org_dir + 'interpro_values.npy', matrix_builder.interpro_values)
            np.save(org_dir + 's2f_values.npy', matrix_builder.s2f_values)
            np.save(org_dir + 'information_content.npy', matrix_builder.information_content)
            goa_values = matrix_builder.goa_values
            hmmer_diff_values = matrix_builder.hmmer_diff_values
            hmmer_cm_values = matrix_builder.hmmer_cm_values
            hmmer_values = matrix_builder.hmmer_values
            interpro_diff_values = matrix_builder.interpro_diff_values
            interpro_cm_values = matrix_builder.interpro_cm_values
            interpro_values = matrix_builder.interpro_values
            s2f_values = matrix_builder.s2f_values
            information_content = matrix_builder.information_content

        if args.compute_metrics:
            # ranges as HX's method
            group_name = ['3-10', '11-30', '31-100', '101-300', '1-300', 'all']
            ranges = [(3, 10), (11, 30), (31, 100), (101, 300), (1, 300), (1, 1000000)]
            group_name = ['all']
            ranges = [(1, 1000000)]

            s2f_nodiff = 0.9 * interpro_values + 0.1 * hmmer_values
            i_h_cm_nodiff = 0.9 * interpro_cm_values + 0.1 * hmmer_cm_values
            sumcol = np.sum(goa_values, axis=0)
            sumrow = np.sum(goa_values, axis=1)

            predictions = []
            gold_standards = []
            information_contents = []
            for low, high in ranges:
                pred = {'I': interpro_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'I+D': interpro_diff_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'I+CM': interpro_cm_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'H': hmmer_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'H+D': hmmer_diff_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'H+CM': hmmer_cm_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'I+H': s2f_nodiff[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'I+H+CM': i_h_cm_nodiff[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                        'I+H+D': s2f_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)]}
                predictions.append(pred)
                gold_standards.append(goa_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)])
                information_contents.append(information_content[(sumcol >= low) & (sumcol <= high)])

            data = defaultdict(list)
            meas = {}
            res = {}
            for pred, gs, r, ic in zip(predictions, gold_standards, group_name, information_contents):
                tell('Processing range:', r)

                # HX measure calculator
                for k in pred.keys():
                    tell(k)
                    rounded_pred = pred[k]
                    if len(np.unique(pred[k])) > 10000:
                        rounded_pred = np.around(pred[k], decimals=4)
                    meas[k] = measures.HX_py(rounded_pred, ic, organism['Tax ID'])
                    tell('computing OVERALL')
                    res[k] = meas[k].compute_overall(gs)
                    for m in res[k].keys():
                        data[m].append([k, r, res[k][m]])

                    tell('computing per-gene')
                    res[k] = meas[k].compute_per_gene(gs)
                    for m in res[k].keys():
                        data[m].append([k, r, res[k][m]])

                    tell('computing per-term')
                    res[k] = meas[k].compute_per_term(gs)
                    for m in res[k].keys():
                        data[m].append([k, r, res[k][m]])

            tell('saving metrics')
            # save all metrics
            metrics = []
            for k, l in data.items():
                for v in l:
                    metrics.append((k, v[0], v[1], v[2]))
            metrics_df = pd.DataFrame(metrics, columns=['metric', 'method', 'range', 'value'])
            metrics_df.to_pickle(org_dir + 'metrics_df.pkl')
    if args.compute_metrics:
        sns.palplot(sns.color_palette(flatui))
        metrics_df = metrics_df.sort_values(by=['metric', 'range', 'method'])
        metrics_df['organism'] = organism['Organism']
        for source, dest in labels:
            metrics_df['method'].replace(source, dest, inplace=True)
        all_organisms = all_organisms.append(metrics_df, ignore_index=True)
        if args.plot_individual:
            tell('Generating figures...')
            for metric in metrics_df['metric'].unique():
                if metric in unsafe_metrics:
                    continue
                plot_selection = metrics_df[metrics_df['metric'] == metric]

                if 'per-' in metric:
                    metric += '(mean)'

                # PLOT WITH RANGES
                a4_dims = (11.7, 8.27)
                fig, ax = plt.subplots(figsize=a4_dims, dpi=100)
                g = sns.barplot(ax=ax, data=plot_selection, x='range', y='value', hue='method', palette=flatui)

                for p in g.patches:
                    ax.annotate("%.3f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                                ha='center', va='center', fontsize=9, xytext=(1, -20), rotation=90,
                                textcoords='offset points', color='white', weight='black')

                ax.set(ylabel=metric)
                ax.set(title=organism['Organism'] + ' ' + organism['Tax ID'])
                ax.legend(bbox_to_anchor=(1.0, 0.5))

                fig.savefig(PLOTS_DIRECTORY+organism['Tax ID']+' '+metric+'-ranges.png')
                plt.close(fig)

                # PLOT WITHOUT RANGES
                fig, ax = plt.subplots(figsize=a4_dims, dpi=100)
                g = sns.barplot(ax=ax, data=plot_selection[plot_selection['range']=='all'],
                                x='method', y='value', palette=flatui)
                for p in g.patches:
                    ax.annotate("%.3f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                                ha='center', va='center', fontsize=11, xytext=(0, 10),
                                textcoords='offset points')
                ax.set(ylabel=metric)
                ax.set(title=organism['Organism'] + ' ' + organism['Tax ID'])
                fig.savefig(PLOTS_DIRECTORY+organism['Tax ID']+' '+metric+'.png')
                plt.close(fig)

if args.compute_metrics:
    tell('Generating figures for all organisms...')
    # plot all organisms
    all_organisms.organism = all_organisms.organism.str.split('(',n=0, expand=True)[0]
    all_organisms.organism = all_organisms.organism.str.split('ATCC',n=0, expand=True)[0]
    all_organisms.organism = all_organisms.organism.str.split('MG',n=0, expand=True)[0]
    all_organisms.to_pickle(os.path.join(MATRICES_DIRECTORY, 'all_organisms.pkl'))
    for metric in all_organisms['metric'].unique():
        if metric in unsafe_metrics:
            continue
        methods = ['InterPro', 'HMMER', 'S2F']
        fig, ax = plt.subplots(figsize=(11.7, 6.27), dpi=100)
        condition = (all_organisms['range'] == 'all') & \
                    (all_organisms['metric'] == metric) & \
                    (all_organisms['method'].isin(methods))
        g = sns.barplot(ax=ax, data=all_organisms[condition],
                        x='organism', y='value', hue='method', palette=flatui3)
        for p in g.patches:
            ax.annotate("%.3f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', fontsize=9, xytext=(1, -20), rotation=90,
                        textcoords='offset points', color='white', weight='black')
        ax.set_xticklabels(ax.get_xticklabels(),rotation=-25, ha='left')
        ax.set(ylabel=metric)
        ax.set(title='')
        ax.legend(bbox_to_anchor=(1.0, 0.5))
        plt.savefig(os.path.join(PLOTS_DIRECTORY, 'ALL '+metric+'.png'), bbox_inches='tight')
        plt.close(fig)
