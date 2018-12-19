"""
Filename: plot_figures.py

Authors:
	Jason Youn - jyoun@ucdavis.edu

Description:

To-do:
"""
import os
import argparse
import numpy as np
import logging as log
import matplotlib.pyplot as plt
from parse_data import ParseData

import time

# default directories
DEFAULT_PARENT_DIR = os.getcwd()
DEFAULT_DATA_DIR = os.path.join(DEFAULT_PARENT_DIR, 'data')
# original data
DEFAULT_PTET_FILENAME = 'pTET.dat'
DEFAULT_PBAD_FILENAME = 'pBAD.dat'
DEFAULT_PBAD_pTET_FILENAME = 'pBAD_pTET.dat'
# predicted data
DEFAULT_PTET_PRED_FILENAME = 'pTET_pred.txt'
DEFAULT_PBAD_DATA_PRED_FILENAME = 'pBAD_pred.txt'
DEFAULT_FINAL_MAT_FILENAME = 'final_mat.txt'

def set_logging():
	"""
	Configure logging.
	"""
	log.basicConfig(format='(%(levelname)s) %(filename)s: %(message)s', level=log.DEBUG)

	# set logging level to WARNING for matplotlib
	logger = log.getLogger('matplotlib')
	logger.setLevel(log.WARNING)

def parse_argument():
	"""
	Parse input arguments.

	Returns:
		- parsed arguments
	"""
	parser = argparse.ArgumentParser(description='Plot figures.')

	parser.add_argument(
		'--data_dir',
		default=DEFAULT_DATA_DIR,
		help='Path to the directory holding all the data.')

	return parser.parse_args()

def plot_actual_with_uncertainty(actual, pred, p_type):
	assert actual.shape == pred.shape
	assert set(list(actual)) == set(list(pred))

	plt.figure(p_type, figsize=(8.5, 7))

	mutant_list = list(actual)[1:]

	if p_type.upper() == 'pTET'.upper():
		x = actual['aTc'].as_matrix()
		plt.xlabel('aTc (ng/ml)')
	elif p_type.upper() == 'pBAD'.upper():
		x = actual['Arabinose'].as_matrix()
		plt.xlabel('Arabinose (%)')
	else:
		raise ValueError('Invalid type!')

	for mutant in mutant_list:
		y_actual = actual[mutant]
		y_pred = pred[mutant]
		error = np.abs(y_actual - y_pred)

		plt.errorbar(x, y_actual, yerr=error, capsize=3)

	plt.xscale('log')
	plt.ylabel('GFP/$OD_{600}$')
	plt.legend(mutant_list)

def plot_double_cascade_comparison(data, case):
	test_cases = list(data)

	# case 1
	case_1 = '10ng/ml aTc'
	case_1_exp = data['{} (exp)'.format(case_1)]
	case_1_sim = data['{} (sim)'.format(case_1)]

	plt.figure(case_1, figsize=(6, 6))
	plt.scatter(case_1_exp, case_1_sim)
	plt.xlabel('{} (exp)'.format(case_1))
	plt.ylabel('{} (sim)'.format(case_1))
	plt.xlim((0,6000))
	plt.ylim((0,6000))

	# case 2
	case_2 = '0.1% ara + 10ng/ml aTc'
	case_2_exp = data['{} (exp)'.format(case_2)]
	case_2_sim = data['{} (sim)'.format(case_2)]

	plt.figure(case_2, figsize=(6, 6))
	plt.scatter(case_2_exp, case_2_sim)
	plt.xlabel('{} (exp)'.format(case_2))
	plt.ylabel('{} (sim)'.format(case_2))
	plt.xlim((0,3000))
	plt.ylim((0,3000))

	# case 3
	case_3 = '100ng/mL aTc'
	case_3_exp = data['{} (exp)'.format(case_3)]
	case_3_sim = data['{} (sim)'.format(case_3)]

	plt.figure(case_3, figsize=(6, 6))
	plt.scatter(case_3_exp, case_3_sim)
	plt.xlabel('{} (exp)'.format(case_3))
	plt.ylabel('{} (sim)'.format(case_3))
	plt.xlim((0,16000))
	plt.ylim((0,16000))


if __name__ == '__main__':
	# set log and parse args
	set_logging()
	args = parse_argument()

	# load original data
	pTET_df = ParseData(os.path.join(args.data_dir, DEFAULT_PTET_FILENAME)).get_dataframe()
	pBAD_df = ParseData(os.path.join(args.data_dir, DEFAULT_PBAD_FILENAME)).get_dataframe()
	pBAD_pTET_df = ParseData(os.path.join(args.data_dir, DEFAULT_PBAD_pTET_FILENAME)).get_dataframe()

	# load predicted data
	pTET_pred_df = ParseData(os.path.join(args.data_dir, DEFAULT_PTET_PRED_FILENAME)).get_dataframe()
	pBAD_pred_df = ParseData(os.path.join(args.data_dir, DEFAULT_PBAD_DATA_PRED_FILENAME)).get_dataframe()
	final_mat_df = ParseData(os.path.join(args.data_dir, DEFAULT_FINAL_MAT_FILENAME)).get_dataframe()

	# plot
	plot_actual_with_uncertainty(pTET_df, pTET_pred_df, 'pTET')
	plot_actual_with_uncertainty(pBAD_df, pBAD_pred_df, 'pBAD')

	plot_double_cascade_comparison(final_mat_df, '10ng/ml aTc')

	plt.show()
