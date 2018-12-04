import pandas as pd
import logging as log

class ParseData(object):
	def __init__(self, filename):
		self.filename = filename
		self._load_data()

	def _load_data(self):
		log.info('Loading file: {}'.format(self.filename))
		self.df = pd.read_csv(self.filename, sep='\t')

	def get_dataframe(self):
		return self.df
