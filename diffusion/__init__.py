import abc

import pandas as pd

from Utils import FancyApp
from seeds import Seed


class Diffusion(FancyApp.FancyApp):

    __metaclass__ = abc.ABCMeta

    @staticmethod
    def _write_results(diffusion_matrix, proteins, terms, filename):
        Seed._seed_to_pandas(diffusion_matrix, 
                             proteins, 
                             terms).to_csv(filename,
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
