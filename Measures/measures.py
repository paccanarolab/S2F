import abc

import numpy as np
import pandas as pd
import sklearn.metrics

from Measures.letor_metrics import ndcg_score
from Utils import FancyApp


class S2FMeasure(FancyApp.FancyApp):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def compute_per_gene(self, real):
        """
        given a matrix with the label samples
        returns the measure result doing it using the 'per gene' mode
        """

    @abc.abstractmethod
    def compute_per_term(self, real):
        """
        given a matrix with the label samples
        returns the measure result doing it using the 'per term' mode
        """

    @abc.abstractmethod
    def compute_overall(self, real):
        """
        given a matrix with the label samples
        returns the measure result doing it using the 'overall' mode
        """

    @abc.abstractmethod
    def setParameters(self, **kwargs):
        '''
        Establishes all the necessary parameters for the measure
        that are defined in showParameters.
        '''

    @abc.abstractmethod
    def showParameters(self):
        """
        Returns the list of parameters that must be passed.
        """


class HX_py(S2FMeasure):

    def __init__(self, prediction, information_content,
                 organism_id, verbose=True):
        super(HX_py, self).__init__()
        self.prediction = prediction
        self.information_content = information_content
        self.__verbose__ = verbose
        self.scores = []
        self.organism_id = organism_id

    def compute_per_gene(self, gold_standard):
        """
        given a matrix with the label samples
        returns the measure result doing it using the 'per gene' mode
        """

        if gold_standard.shape[0] < 1:
            return {}
        # we iterate all the rows and compute the mean of the row-wise scores
        genewise = []
        for row in range(gold_standard.shape[0]):
            hx = HX_py.HX_iteration(
                self.prediction[row, :],
                gold_standard[row, :]
            )
            th_list, idx = np.unique(self.prediction[row, :],
                                     return_index=True)
            s = self.compute_s_measure(self.prediction[row, :],
                                       gold_standard[row, :], th_list)
            genewise.append({**hx, **s})
        # rename all the keys
        keys = genewise[0].keys()
        result = {}
        valid_keys = [k for k in keys if k not in
                      ['s', 'ru', 'mi', 'roc', 'pr']]
        for k in valid_keys:
            result[k+' per-gene'] = np.mean([g[k] for g in genewise])
            result[k+' per-gene raw'] = np.array([g[k] for g in genewise])
        return result

    def compute_per_term(self, gold_standard):
        """
        given a matrix with the label samples
        returns the measure result doing it using the 'per term' mode
        """
        if gold_standard.shape[1] < 1:
            return {}
        # we iterate all the columns and compute the mean
        # of the column-wise scores
        termwise = []
        for col in range(gold_standard.shape[1]):
            hx = HX_py.HX_iteration(
                self.prediction[:, col],
                gold_standard[:, col]
            )
            th_list, idx = np.unique(self.prediction[:, col],
                                     return_index=True)
            s = self.compute_s_measure(
                self.prediction[:, col], gold_standard[:, col], th_list,
                ic=self.information_content[col])
            termwise.append({**hx, **s})
        # rename all the keys
        keys = termwise[0].keys()
        result = {}
        valid_keys = [k for k in keys if k not in
                      ['s', 'ru', 'mi', 'roc', 'pr']]
        for k in valid_keys:
            result[k + ' per-term'] = np.mean([g[k] for g in termwise])
            result[k + ' per-term raw'] = np.array([g[k] for g in termwise])
        return result

    def compute_overall(self, gold_standard):
        """
        given a matrix with the label samples
        returns the measure result doing it using the 'overall' mode
        """
        hx = HX_py.HX_iteration(
            self.prediction.flatten(),
            gold_standard.flatten(),
        )

        th_list, idx = np.unique(self.prediction.flatten(), return_index=True)
        s = self.compute_s_measure(self.prediction, gold_standard, th_list)
        return {**hx, **s}

    def compute_s_measure(self, pred, gs, ths, ic='infer'):
        # masked = gs * pred
        ru = []
        mi = []
        s = []
        smin = np.inf
        for th in ths:
            false_negatives = ((pred < th) & (gs > 0)).astype(np.float)
            false_positives = ((pred >= th) & (gs < 1)).astype(np.float)
            if type(ic) == str and ic == 'infer':
                # TODO: update this to positive
                ic = -self.information_content
            ru.append(np.sum(false_negatives * ic)/gs.shape[0])
            mi.append(np.sum(false_positives * ic)/gs.shape[0])
            s.append(np.sqrt(ru[-1]**2 + mi[-1]**2))
            if s[-1] < smin:
                smin = s[-1]
        return {
            's': s,
            'ru': ru,
            'mi': mi,
            'smin': smin,
        }

    @staticmethod
    def HX_iteration(prediction, gold_standard, recall_th=0.2):
        N = np.prod(gold_standard.shape)
        dat = np.zeros((N, 2))
        dat[:, 0] = prediction.flatten()
        dat[:, 1] = gold_standard.flatten()
        dat = dat[dat[:, 0].argsort()]

        qty_pos = np.sum(dat[:, 1] > 0)
        qty_neg = N - qty_pos

        # threshold list
        th_list, idx = np.unique(dat[:, 0], return_index=True)

        qty_pos = np.sum(dat[:, 1] > 0)
        qty_neg = N - qty_pos

        n_steps = len(th_list)

        # if it does not make any sense to compute all the score,
        # return HX's 'default' values
        if qty_pos * qty_neg == 0 or n_steps < 2:
            auc = 0.5
            pre_at = 0
            # f1 = 2 / (1 + len(dat[:, 1]) /
            #           (qty_pos + np.finfo(np.double).tiny))
            ndcg = 0
            jac = 0
            if qty_pos == 0:
                jac = 0
            if qty_neg == 0:
                jac = 1
            if n_steps > 0:
                jac = qty_pos / len(dat[:, 1])
            return {
                'AUC': auc,
                # based on https://stats.stackexchange.com/a/266989/227894
                'AUPR': qty_pos / N,
                'Precision at 0.2 Recall': pre_at,
                'F_max': 0,
                'NDCG': ndcg,
                'Jaccard': jac,
                'roc': [np.array([0, 1]), np.array([0, 1])],
                'pr': [np.array([0, 1]), np.array([qty_pos / N, qty_pos / N])],
            }

        FP, TP, TH = sklearn.metrics.ranking._binary_clf_curve(dat[:, 1],
                                                               dat[:, 0],
                                                               pos_label=1.0)

        # TN = qty_neg - FP
        FN = qty_pos - TP

        FPR = FP / qty_neg
        TPR = TP / qty_pos

        rec = TP / (TP + FN)
        pre = TP / (TP + FP)

        f_max = np.max(2.0 * (np.multiply(pre, rec)) /
                             (pre + rec + np.finfo(np.double).tiny))

        auc = sklearn.metrics.auc(FPR, TPR)
        aupr = sklearn.metrics.auc(rec, pre)

        jac = np.max(TP / (TP + FN + FP))

        pre_at_idx = np.where(rec >= recall_th)
        pre_at_idx = np.argmin(np.abs(rec[pre_at_idx] - recall_th))
        pre_at = pre[pre_at_idx]

        ndcg = ndcg_score(gold_standard, prediction, k=N)
        # vec = np.log2(np.arange(1, N + 1))
        # vec[0] = 1

        return {
            'AUC': auc,
            'AUPR': aupr,
            'Precision at 0.2 Recall': pre_at,
            'F_max': f_max,
            'NDCG': ndcg,
            'Jaccard': jac,
            'roc': [FPR, TPR],
            'pr': [rec, pre]
            # 'tp': TP,
            # 'fp': FP,
            # 'tn': TN,
            # 'fn': FN,
        }

    def setParameters(self, **kwargs):
        """
        Establishes all the necessary parameters for the measure
        that are defined in showParameters.
        """
        try:
            self.prediction = kwargs['prediction']
        except KeyError:
            pass

        try:
            self.evidence_codes = kwargs['evidence_codes']
        except KeyError:
            pass

    def showParameters(self):
        """
        prediction: numpy array that contains the predictions made by the
            PFP method This code assumes that the shape of the matrix has
            the genes in the rows and the terms in the columns
        evidence_codes: to filter the GOA by evidence codes, this is the
            list of evidence codes to be considered as true
        """
        return ['prediction', 'evidence_codes']


def log(api, message):
    try:
        api.send_direct_message(screen_name='torresmateo', text=message)
    except Exception:
        print('could not tweet')
    print(message)


if __name__ == '__main__':

    import seaborn as sns
    import matplotlib.pyplot as plt
    from Utils.notification import load_configuration, get_api
    import os.path
    import _pickle as pickle
    home = os.path.expanduser('~')
    api = get_api(load_configuration(home+'/.config/pcva036notification.json'))

    organisms = [
        '1111708',
        '1392',
        '208964',
        '223283',
        '224308',
        '243277',
        '246200',
        '83332',
        '99287',
        '83333'
    ]

    plot_data = {}

    for org in organisms:

        log(api, '[S2F] generating plot information for organism ' + org)
        goa_values = np.load('matrices/' + org + '_goa_values.npy')
        hmmer_diff_values = np.load('matrices/' + org +
                                    '_hmmer_diff_values.npy')
        hmmer_values = np.load('matrices/' + org + '_hmmer_values.npy')
        interpro_diff_values = np.load('matrices/' + org +
                                       '_interpro_diff_values.npy')
        interpro_values = np.load('matrices/' + org + '_interpro_values.npy')
        s2f_values = np.load('matrices/' + org + '_s2f_values.npy')

        group_name = ['3-10', '11-30', '31-100', '101-300', '1-300']

        s2f_nodiff = 0.9*interpro_values + 0.1*hmmer_values
        sumcol = np.sum(goa_values, axis=0)
        sumrow = np.sum(goa_values, axis=1)

        ranges = [(3, 10), (11, 30), (31, 100), (101, 300), (1, 300)]
        predictions = []
        gold_standards = []
        for low, high in ranges:
            pred = {}
            pred['I'] = interpro_values[sumrow > 3, :][:, (sumcol >= low) &
                                                       (sumcol <= high)]
            pred['I+D'] = interpro_diff_values[sumrow > 3, :][:,
                                                              (sumcol >= low) &
                                                              (sumcol <= high)]
            pred['H'] = hmmer_values[sumrow > 3, :][:, (sumcol >= low) &
                                                       (sumcol <= high)]
            pred['H+D'] = hmmer_diff_values[sumrow > 3, :][:, (sumcol >= low) &
                                                              (sumcol <= high)]
            pred['I+H'] = s2f_nodiff[sumrow > 3, :][:, (sumcol >= low) &
                                                       (sumcol <= high)]
            pred['I+H+D'] = s2f_values[sumrow > 3, :][:, (sumcol >= low) &
                                                         (sumcol <= high)]
            predictions.append(pred)
            gold_standards.append(goa_values[sumrow > 3, :][:, (sumcol >= low)
                                                            & (sumcol <= high)]
                                  )

        data = {
            '1 - AUC Overall': [], '1 - AUC Per Gene': [],
            '1 - AUC Per Term': [], 'F1 Score Overall': [],
            'F1 Score Per Gene': [], 'F1 Score Per Term': [],
        }

        for pred, gs, r in zip(predictions, gold_standards, group_name):
            print('Processing range:', r)
            # HX measure calculator
            meas_i = HX_py(pred['I'], org + '_I_' + r)
            meas_id = HX_py(pred['I+D'], org + '_I+D_' + r)
            meas_h = HX_py(pred['H'], org + '_H_')
            meas_hd = HX_py(pred['H+D'], org + '_H+D_' + r)
            meas_ih = HX_py(pred['I+H'], org + '_I+H_' + r)
            meas_ihd = HX_py(pred['I+H+D'], org + '_I+H+D_' + r)

            log(api, org + ' InterPro...' + r)
            res_i = meas_i.compute_overall(gs)
            data['1 - AUC Overall'].append(['I', r, 1 - res_i['auc']])

            log(api, org + ' InterPro + Diffusion...' + r)
            res_id = meas_id.compute_overall(gs)
            data['1 - AUC Overall'].append(['I+D', r, 1 - res_id['auc']])

            log(api, org + ' HMMer...' + r)
            res_h = meas_h.compute_overall(gs)
            data['1 - AUC Overall'].append(['H', r, 1 - res_h['auc']])

            log(api, org + ' HMMer + Diffusion...' + r)
            res_hd = meas_hd.compute_overall(gs)
            data['1 - AUC Overall'].append(['H+D', r, 1 - res_hd['auc']])

            log(api, org + ' InterPro + HMMer...' + r)
            res_ih = meas_ih.compute_overall(gs)
            data['1 - AUC Overall'].append(['I+H', r, 1 - res_ih['auc']])

            log(api, org + ' S2F...' + r)
            res_ihd = meas_ihd.compute_overall(gs)
            data['1 - AUC Overall'].append(['I+H+D', r, 1 - res_ihd['auc']])

            plot_data[org+r] = {'I': res_i,
                                'I+D': res_id,
                                'H': res_h,
                                'H+D': res_hd,
                                'I+H': res_ih,
                                'I+H+D': res_ihd}

        print(data)
        auc_overall = pd.DataFrame(data['1 - AUC Overall'],
                                   columns=['method', 'range', 'auc'])

        a4_dims = (11.7, 8.27)
        fig, ax = plt.subplots(figsize=a4_dims)
        sns.barplot(ax=ax, data=auc_overall, x='range', y='auc', hue='method')
        ax.set(ylabel='1 - AUC')
        plt.savefig('1-AUC' + org + '.png')
        log(api, '[S2F] organism ' + org + ' done')
    log(api, '[S2F] plot script done')

    with open('plot_data.p', 'wb') as fp:
        pickle.dump(plot_data, fp)
