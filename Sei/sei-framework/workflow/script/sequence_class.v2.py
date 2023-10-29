"""
Description:
    This script is run after `vep_cli.py`. It computes the sequence-class
    level variant effect scores based on Sei chromatin profile predictions,
    sorts the variants based on the maximum absolute scores across sequence classes
    and outputs the results as TSV files.

Usage:
    sequence_class.py <results-dir> <vcf>
    sequence_class.py -h | --help

Options:
    <results-dir>    The results directory.
    <vcf>            Name of VCF file.

"""
import os

from docopt import docopt
import h5py
import numpy as np
import pandas as pd


def get_targets(filename):
    targets = []
    with open(filename, 'r') as file_handle:
        for line in file_handle:
            targets.append(line.strip())
    return targets


def get_data(filename):
    fh = h5py.File(filename, 'r')
    data = fh["data"][()]
    fh.close()
    return data

if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    results_dir = arguments['<results-dir>']
    input_vcf = arguments['<vcf>']

    _, filename = os.path.split(input_vcf)
    filename_prefix = '.'.join(filename.split('.')[:-1])

    sei_dir = "/home/sonic/sei/sei-framework/model"

    profile_pred_dir = os.path.join(results_dir, 'chromatin-profiles-hdf5')

    chromatin_profile_ref = get_data(os.path.join(
        profile_pred_dir,
        "{0}.ref_predictions.h5".format(filename_prefix)))
    chromatin_profile_alt = get_data(os.path.join(
        profile_pred_dir,
        "{0}.alt_predictions.h5".format(filename_prefix)))

    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_targets.npy'))
    histone_inds = np.load(os.path.join(sei_dir, 'histone_inds.npy'))

    chromatin_profile_ref_adjust = chromatin_profile_ref.copy()
    chromatin_profile_ref_adjust[:, histone_inds] = \
        chromatin_profile_ref_adjust[:, histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_ref[:, histone_inds], axis=1))[:, None]

    chromatin_profile_alt_adjust = chromatin_profile_alt.copy()
    chromatin_profile_alt_adjust[:, histone_inds] = \
        chromatin_profile_alt_adjust[:, histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_alt[:, histone_inds], axis=1))[:, None]

    refproj = (np.dot(chromatin_profile_ref_adjust,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))
    altproj = (np.dot(chromatin_profile_alt_adjust,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))

    np.savetxt(os.path.join(results_dir, 'sequence_class_ref_scores.tsv.gz'), refproj, delimiter='\t')
    np.savetxt(os.path.join(results_dir, 'sequence_class_alt_scores.tsv.gz'), altproj, delimiter='\t')
    
