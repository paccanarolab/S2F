from Measures import BuildMatrices, measures
import numpy as np
import pandas as pd
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
from Utils.notification import load_configuration, get_api

dataroot = '/home/paccanaro/COMMON/PROJECTS/S2F/data/2017/'
obo = dataroot + 'INPUT/go.obo'
DATA_DIRECTORY =  '../../../data/bacteria_selection_2017/'
GOA_DIRECTORY = DATA_DIRECTORY + 'selection_goa/'
FASTA_DIRECTORY = DATA_DIRECTORY + 'selected_fasta/'
PLOTS_DIRECTORY = DATA_DIRECTORY + 'comparison_plots/'
MATRICES_DIRECTORY = PLOTS_DIRECTORY + 'matrices/'

#TODO: once the S2F run is finished, this list will be useless
finished_list = [
    '1111708',
    '1392',
    '208964',
    '223283',
    '224308',
    '243277',
    '246200',
    '83332',
    '83333',
    '99287']

for organism in finished_list:
    interpro_file = dataroot + 'seeds/interpro/'+organism+'.matrix'
    hmmer_file = dataroot + 'seeds/hmmer/'+organism+'.matrix'
    for f in os.listdir(dataroot + 'output/'):
        if f.startswith(organism):
            interpro_diff_file = dataroot + 'output/' + f + '/interpro.diffusion'
            hmmer_diff_file = dataroot + 'output/' + f + '/hmmer.diffusion'
            s2f_file = dataroot + 'output/' + f + '/final.diffusion'
    goa = dataroot + 'INPUT/goa/' + organism + '.goa'
    fasta = dataroot + 'INPUT/sequences/' + organism + '.fasta'
    obo = dataroot + 'INPUT/go.obo'

    # load or build matrices
    org_dir = MATRICES_DIRECTORY + organism + '/'
    if os.path.exists(org_dir):
        goa_values = np.load(org_dir + 'goa_values.npy')
        hmmer_diff_values = np.load(org_dir + 'hmmer_diff_values.npy')
        hmmer_values = np.load(org_dir + 'hmmer_values.npy')
        interpro_diff_values = np.load(org_dir + 'interpro_diff_values.npy')
        interpro_values = np.load(org_dir + 'interpro_values.npy')
        s2f_values = np.load(org_dir + 's2f_values.npy')
    else:
        os.makedirs(org_dir)
        matrix_builder = BuildMatrices.BuildMatrices(
            interpro_file, interpro_diff_file, hmmer_file, hmmer_diff_file, s2f_file, goa, fasta, obo)
        np.save(org_dir + 'goa_values.npy', matrix_builder.goa_values)
        np.save(org_dir + 'hmmer_diff_values.npy', matrix_builder.hmmer_diff_values)
        np.save(org_dir + 'hmmer_values.npy', matrix_builder.hmmer_values)
        np.save(org_dir + 'interpro_diff_values.npy', matrix_builder.interpro_diff_values)
        np.save(org_dir + 'interpro_values.npy', matrix_builder.interpro_values)
        np.save(org_dir + 's2f_values.npy', matrix_builder.s2f_values)
        goa_values = matrix_builder.goa_values
        hmmer_diff_values = matrix_builder.hmmer_diff_values
        hmmer_values = matrix_builder.hmmer_values
        interpro_diff_values = matrix_builder.interpro_diff_values
        interpro_values = matrix_builder.interpro_values
        s2f_values = matrix_builder.s2f_values

    # ranges as HX's method
    group_name = ['3-10', '11-30', '31-100', '101-300', '1-300', 'all']

    s2f_nodiff = 0.9 * interpro_values + 0.1 * hmmer_values
    sumcol = np.sum(goa_values, axis=0)
    sumrow = np.sum(goa_values, axis=1)

    ranges = [(3, 10), (11, 30), (31, 100), (101, 300), (1, 300), (1, 1000000)]
    predictions = []
    gold_standards = []
    for low, high in ranges:
        pred = {'I': interpro_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                'I+D': interpro_diff_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                'H': hmmer_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                'H+D': hmmer_diff_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                'I+H': s2f_nodiff[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)],
                'I+H+D': s2f_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)]}
        predictions.append(pred)
        gold_standards.append(goa_values[sumrow > 3, :][:, (sumcol >= low) & (sumcol <= high)])

    data = {
        '1 - AUC Overall': [], '1 - AUC Per Gene': [], '1 - AUC Per Term': [],
        'F1 Score Overall': [], 'F1 Score Per Gene': [], 'F1 Score Per Term': [],

    }
    if os.path.exists(org_dir + 'metrics_df.pkl'):
        metrics_df = pd.read_pickle(org_dir + 'metrics_df.pkl')
    else:
        data = {
            'AUC Overall': [], 'AUC Overall sklearn': [], 'AUC Per Gene': [], 'AUC Per Term': [],
            'F1 Score Overall': [], 'F1 Score Per Gene': [], 'F1 Score Per Term': [],

        }
        for pred, gs, r in zip(predictions, gold_standards, group_name):
            print('Processing range:', r)
            # HX measure calculator
            meas_i = measures.HX_py(pred['I'], organism)
            meas_id = measures.HX_py(pred['I+D'], organism)
            meas_h = measures.HX_py(pred['H'], organism)
            meas_hd = measures.HX_py(pred['H+D'], organism)
            meas_ih = measures.HX_py(pred['I+H'], organism)
            meas_ihd = measures.HX_py(pred['I+H+D'], organism)

            data['AUC Overall'].append(['I', r, meas_i.compute_overall(gs)['auc']])
            data['AUC Overall'].append(['I+D', r, meas_id.compute_overall(gs)['auc']])
            data['AUC Overall'].append(['H', r, meas_h.compute_overall(gs)['auc']])
            data['AUC Overall'].append(['H+D', r, meas_hd.compute_overall(gs)['auc']])
            data['AUC Overall'].append(['I+H', r, meas_ih.compute_overall(gs)['auc']])
            data['AUC Overall'].append(['I+H+D', r, meas_ihd.compute_overall(gs)['auc']])

            meas_i = measures.OneAUC(pred['I'], organism)
            meas_id = measures.OneAUC(pred['I+D'], organism)
            meas_h = measures.OneAUC(pred['H'], organism)
            meas_hd = measures.OneAUC(pred['H+D'], organism)
            meas_ih = measures.OneAUC(pred['I+H'], organism)
            meas_ihd = measures.OneAUC(pred['I+H+D'], organism)

            data['AUC Overall sklearn'].append(['I', r, meas_i.compute_overall(gs)])
            data['AUC Overall sklearn'].append(['I+D', r, meas_id.compute_overall(gs)])
            data['AUC Overall sklearn'].append(['H', r, meas_h.compute_overall(gs)])
            data['AUC Overall sklearn'].append(['H+D', r, meas_hd.compute_overall(gs)])
            data['AUC Overall sklearn'].append(['I+H', r, meas_ih.compute_overall(gs)])
            data['AUC Overall sklearn'].append(['I+H+D', r, meas_ihd.compute_overall(gs)])

        # save all metrics
        metrics = []
        for k, l in data.items():
            for v in l:
                metrics.append((k, v[0], v[1], v[2]))
        metrics_df = pd.DataFrame(metrics, columns=['metric', 'method', 'range', 'value'])
        metrics_df.to_pickle(org_dir + 'metrics_df.pkl')

    flatui = ["#705319", "#e5aa34", "#911f76", "#e534bb", "#1c7c69", "#31d6b5", ]
    sns.palplot(sns.color_palette(flatui))
    for metric in metrics_df['metric'].unique():
        plot_selection = metrics_df[metrics_df['metric'] == metric]

        a4_dims = (11.7, 8.27)
        fig, ax = plt.subplots(figsize=a4_dims, dpi=100)
        g = sns.barplot(ax=ax, data=plot_selection, x='range', y='value', hue='method', palette=flatui)
        ax.set(ylabel='AUC')
        ax.set(title=organism)
        fig.savefig(PLOTS_DIRECTORY+organism+' '+metric+'-ranges.png')

        fig, ax = plt.subplots(figsize=a4_dims, dpi=100)
        g = sns.barplot(ax=ax, data=plot_selection[plot_selection['range'] == 'all'],
                        x='method', y='value', palette=flatui)
        for p in g.patches:
            ax.annotate("%.3f" % p.get_height(), (p.get_x() + p.get_width() / 2., p.get_height() - 0.02),
                        ha='center', va='center', fontsize=11, xytext=(0, 20),
                        textcoords='offset points')
        ax.set(ylabel='AUC')
        ax.set(title=organism)
        fig.savefig(PLOTS_DIRECTORY+organism+' '+metric+'.png')
        plt.close(fig)