import sys

import numpy as np
import pandas as pd
import scipy
import subprocess
import os


def all_indices(value, qlist):
    indices = []
    idx = -1
    while True:
        try:
            idx = qlist.index(value, idx+1)
            indices.append(idx)
        except ValueError:
            break

    return indices


# taken from https://stackoverflow.com/questions/3041986/\
    # apt-command-line-interface-like-yes-no-input
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def jaccard(A):
    """ Makes a similarity matrix NxN from
        an NxM matrix, by multiplying it with itself
        and dividing by a normalising constant
        :param A: scipy sparse matrix of shape (N, M)
    """
    N, M = A.shape
    if N == M:
        # I don't understand this and likely,
        # will never be executed
        A = A - np.diag(np.diag(A)) + np.eye(N)

    Co = A.dot(A.transpose())
    CLEN = np.tile(A.sum(axis=1).ravel(), (N, 1))

    return Co / (CLEN + CLEN.transpose() - Co + 1)


# https://gist.github.com/bwhite/3726239#gistcomment-1509098
def dcg_at_k(r, k):
    r = np.asfarray(r)[:k]
    if r.size:
        return np.sum(np.subtract(np.power(2, r), 1) /
                      np.log2(np.arange(2, r.size + 2)))
    return 0.


def ndcg_at_k(r, k):
    idcg = dcg_at_k(sorted(r, reverse=True), k)
    if not idcg:
        return 0.
    return dcg_at_k(r, k) / idcg


def extract_uniprot_accession(protein_id):
    try:
        return protein_id.split('|')[1]
    except IndexError:
        return protein_id


def keep_entire_prot_id(fasta_line):
    return fasta_line[1:].split()[0]


def keep_uniprot_accession(fasta_line):
    try:
        return extract_uniprot_accession(fasta_line[1:].split()[0])
    except IndexError:
        return keep_entire_prot_id(fasta_line)


def extract_indices_from_fasta(fasta, processing_func=keep_entire_prot_id):
    """
    given a fasta file, it will build a dataframe with
    the proteins as index, very useful for indexing graphs
    :param fasta: fasta file
    :param processing_func: the function to process the line in case
            the entire ID is not required
    :return: pandas.DataFrame
    """
    proteins = []
    for line in open(fasta):
        if line.startswith('>'):
            proteins.append(processing_func(line))
    proteins_df = pd.DataFrame(list(enumerate(sorted(proteins))))
    proteins_df.columns = ['protein idx', 'protein id']
    proteins_df.set_index('protein id', inplace=True)
    return proteins_df
    

def line_count(filename):
    if os.name == 'nt':  # windows
        # apparently the best way to count lines in
        # windows: https://superuser.com/a/959037/1105734
        out = subprocess.check_output(r'find /c /v "" {f}'.format(f=filename))
        return int(out.split(b' ')[-1])
    else:  # unix
        # taken from
        # https://gist.github.com/zed/0ac760859e614cd03652#file-gistfile1-py-L41
        out = subprocess.Popen(['wc', '-l', filename],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT).communicate()[0]
        return int(out.partition(b' ')[0])
