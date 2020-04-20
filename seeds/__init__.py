import abc
import pandas as pd
from Utils import FancyApp


class Seed(FancyApp.FancyApp):

    __metaclass__ = abc.ABCMeta

    @staticmethod
    def _seed_to_pandas(seed_matrix, proteins, terms):
        """
        transform a seed matrix into a pandas DataFrame
        Parameters
        ----------
        seed_matrix : scipy.sparse matrix
            The labelling
        proteins : pandas.DataFrame
            Used to map the indices on the `seed_matrix` to the protein ids
        terms:
            Used to map the indices on the `seed_matrix` to the GO ids
        Returns
        -------
        pandas.DataFrame
        """
        data = {
            'protein': seed_matrix.row,
            'goterm': seed_matrix.col,
            'score': seed_matrix.data
        }
        labelling_df = pd.DataFrame(data)
        labelling_df = labelling_df.merge(proteins.reset_index(),
                                          left_on='protein',
                                          right_on='protein idx')
        labelling_df = labelling_df.merge(terms.reset_index(),
                                          left_on='goterm',
                                          right_on='term idx')
        return labelling_df[['protein id', 'term id', 'score']]

    @classmethod
    def get_seed(cls):
        """
        Returns a matrix with the seed ready to be diffused
        Returns
        -------
        scipy.sparse matrix
            a matrix with shape `(n, m)` where `n` is the number of
            proteins and `m` the number of terms
        """

    @classmethod
    def process_output(cls, **kwargs):
        """
        Parses the output of the prediction method and extracts its annotations
        Parameters
        ----------
        kwargs:
            To be handled by the particular needs of the implementation.
        """
