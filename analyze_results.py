"""
Filename: analyze_results.py

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
from pathlib import Path
from parse_data import ParseData
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

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
# parameter tuning result
DEFAULT_PARAM_TUNING_FILENAME = 'param_tuning_results.txt'

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
	parser = argparse.ArgumentParser(description='Analyze the simulation results.')

	parser.add_argument(
		'--data_dir',
		metavar='data directory',
		default=DEFAULT_DATA_DIR,
		help='Path to the directory holding all the data.')

	parser.add_argument(
		'--show_plots',
		action='store_true',
		default=False,
		help='Display plots if set.')

	parser.add_argument(
		'--param_search',
		action='store_true',
		default=False,
		help='Run in parameter search mode if set.')

	parser.add_argument(
		'-k',
		'--k_convert_gfp_to_tetR',
		metavar='conversion constant',
		type=int,
		help='Conversion constant to convert GFP to tetR')

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

def plot_double_cascade_comparison(data):
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

def plotROC(exp, sim, show_plots):
	fpr, tpr , _ = roc_curve(exp, sim)
	auROC = auc(fpr, tpr)

	if show_plots:
		# plot figure
		fig = plt.figure('ROC')
		ax1 = fig.add_subplot(111)
		ax1.plot(fpr, tpr, lw=2, color='darkorange', label='AUC:{:.3f}'.format(auROC))
		ax1.plot([0, 1], [0, 1], color='navy', linestyle='--')
		plt.title('ROC')
		plt.xlabel('False Positive Rate')
		plt.ylabel('True Positive Rate')
		plt.legend(loc='lower right')
		# plt.savefig('fig/ROC.pdf')

	return auROC

def plotPR(exp, sim, show_plots):
	precision, recall , _ = precision_recall_curve(exp.ravel(), sim.ravel())
	avg_prec = average_precision_score(exp, sim, average='micro')

	if show_plots:
		# plot figure
		fig = plt.figure('PR')
		ax1 = fig.add_subplot(111)
		ax1.plot(recall, precision, lw=2, color='darkorange', label='AP:{:.3f}'.format(avg_prec))
		plt.title('PR')
		plt.xlabel('Recall')
		plt.ylabel('Precision')
		plt.legend(loc='upper right')
		# plt.savefig('PR.pdf')

	return avg_prec

def save_param_search(data_dir, k_convert_gfp_to_tetR, threshold, num_tunable_exp, auroc):
	file_path = Path(os.path.join(data_dir, DEFAULT_PARAM_TUNING_FILENAME))

	if file_path.exists() is False:
		with open(file_path, 'w+') as file:
			file.write('k_convert_gfp_to_tetR\tthreshold\tnum_tunable_exp\tauroc\n')

	with open(file_path, 'a') as file:
		file.write('{}\t{}\t{}\t{}\n'.format(
			k_convert_gfp_to_tetR,
			threshold,
			num_tunable_exp,
			auroc))

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
	if args.show_plots:
		plot_actual_with_uncertainty(pTET_df, pTET_pred_df, 'pTET')
		plot_actual_with_uncertainty(pBAD_df, pBAD_pred_df, 'pBAD')
		plot_double_cascade_comparison(final_mat_df)

	# find folds to be used for label for (exp)
	ara_fold_exp = final_mat_df.loc[:, '10ng/ml aTc (exp)'] / final_mat_df.loc[:, '0.1% ara + 10ng/ml aTc (exp)']

	# find folds to be used for label for (sim)
	ara_fold_sim = final_mat_df.loc[:, '10ng/ml aTc (sim)'] / final_mat_df.loc[:, '0.1% ara + 10ng/ml aTc (sim)']

	thres_linspace = np.linspace(ara_fold_exp.min(), ara_fold_exp.max(), num=20)[1:-1]

	for threshold in thres_linspace:
		label_exp = np.where(ara_fold_exp > threshold, 'Tunable', 'non-tunable')

		unique_exp, unique_counts = np.unique(label_exp, return_counts=True)
		if unique_exp.shape[0] != 2:
			continue

		num_tunable_exp = int(unique_counts[np.where(unique_exp == 'Tunable')])

		# convert labels to one-hot vectors
		label_exp_one_hot = label_binarize(label_exp, unique_exp).ravel()

		# analysis
		# avg_prec = plotPR(label_exp_one_hot, label_sim_one_hot)
		auroc = plotROC(label_exp_one_hot, ara_fold_sim.values, args.show_plots)

		save_param_search(args.data_dir, args.k_convert_gfp_to_tetR, threshold, num_tunable_exp, auroc)

	# only show plots if requested
	if args.show_plots: plt.show()
