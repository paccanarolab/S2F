import abc

from Utils import FancyApp


class Seed(FancyApp.FancyApp):

    __metaclass__ = abc.ABCMeta

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
