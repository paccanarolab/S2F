import abc

import pandas as pd

from Utils import FancyApp


class Diffusion(FancyApp.FancyApp):

    __metaclass__ = abc.ABCMeta

    @staticmethod
    def _write_results(diffusion_matrix, proteins, terms, filename):
        data = {
            'protein': diffusion_matrix.row,
            'goterm': diffusion_matrix.col,
            'score': diffusion_matrix.data
        }
        labelling_df = pd.DataFrame(data)
        labelling_df = labelling_df.merge(proteins.reset_index(),
                                          left_on='protein',
                                          right_on='protein idx')
        labelling_df = labelling_df.merge(terms.reset_index(),
                                          left_on='goterm',
                                          right_on='term idx')
        labelling_df[['protein id', 'term id', 'score']].to_csv(filename,
                                                                sep='\t',
                                                                index=False,
                                                                header=None)

    @abc.abstractmethod
    def write_results(self, filename):
        """
        Writes the diffusion results into a file in TSV format
        :param filename: path to the desired output file
        """

    @abc.abstractmethod
    def diffuse(self, initial_guess, **kwargs):
        """
        Diffuses the labels found in `initial_guess` into the computed kernel

        Parameters
        ----------
        initial_guess : scipy.sparse matrix
            The labelling that should be propagated in the graph
        kwargs
            will be handled by the implementation accordingly

        Returns
        -------
        scipy.sparse matrix
            The final labelling after propagation
        """

    @abc.abstractmethod
    def compute_kernel(self, **kwargs):
        """
        If the diffusion method should make any transformations on the graph,
        it should be done in this function

        Parameters
        ----------
        kwargs
            will be handled by the implementation accordingly
        """
