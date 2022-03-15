import os
import sys
import numpy as np
import pandas as pd
from tools.io import *
from tdc_metrics import calculate_correlation_make_plots
from sklearn import linear_model
from sklearn.linear_model import ElasticNet
from sklearn.exceptions import ConvergenceWarning
import scipy.stats

import warnings
warnings.simplefilter("ignore", ConvergenceWarning)
warnings.simplefilter("ignore", RuntimeWarning)
print('Ignoring any warnings of type: ConvergenceWarning, RuntimeWarning')

class TumorDeconvolutionModels:

    def __init__(self, model_level, model_name = 'arithmetic_mean', signature_matrix = None, round = None,
                 input_name = None, df_data = None, gs = None, use_gold_standard = True, recordCollectionToBinaries = False):

        self.recordCollectionToBinaries = recordCollectionToBinaries

        self.input_name = input_name
        self.model_level = model_level
        self.model_name = model_name
        self.round = round
        self.inputDir = os.path.join('input', '')
        self.outputDir = os.path.join('output', '')
        self.df_gold_standard = gs
        self.df_data = df_data

        #if signature_matrix is None:
        #    signature_matrix = None

        self.df_signature = pd.read_excel(signature_matrix, index_col=0,
                                     header=[0,1] if model_level == 'fine' else [0],
                                     skiprows=[2] if model_level == 'fine' else [1])

        if self.df_data is None:
            if self.input_name is None:
                print('Provide input file name')
                return
            else:
                self.df_data = self.readDatasets()

        if use_gold_standard:
            if self.df_gold_standard is None:
                if self.round is None:
                    print('Specify round')
                    return
                else:
                    try:
                        self.df_gold_standard = pd.read_csv('input/gold_standards/lb_%s_r%s.csv' % (model_level, round),
                                              index_col=[0,1,2], header=0)['measured']
                    except Exception as exception:
                        print(exception)

        return

    def readDatasets(self, printDetail = False):

        # Create and populate pd.DataFrame 'df_data'.
        # Rows are genes, columns are sample and dataset names
        df_meta_data = pd.read_csv(os.path.join(self.inputDir, self.input_name))

        df_data = pd.DataFrame()

        for row, s_meta_data in df_meta_data.iterrows():
            print('Loading dataset:', s_meta_data.loc['dataset.name'])

            if printDetail:
                for key, value in s_meta_data.iloc[1:].iteritems():
                    print('\t', key, ': ', value)

            # Get dataset name
            dataset = s_meta_data.loc['dataset.name']

            # Read data
            df_data_temp = pd.read_csv(os.path.join(self.inputDir, s_meta_data['hugo.expr.file'].strip('"').strip('/').strip('input').strip('/')), index_col='Gene')

            # Replace missing values with 0
            df_data_temp.fillna(0, inplace=True)

            # Convert data to linear scale
            if s_meta_data.loc['scale'] == 'linear':
                print('Data scale is linear')
            elif s_meta_data.loc['scale'] == 'log2':
                print('Data scale is log2. Converting to linear')
                df_data_temp = np.power(2., df_data_temp) - 1.
            elif s_meta_data.loc['scale'] == 'ln':
                print('Data scale is ln.  Converting to linear')
                df_data_temp = np.exp(df_data_temp)
            elif s_meta_data.loc['scale'] == 'log10':
                print('Data scale is log10. Converting to linear')
                df_data_temp = np.power(10., df_data_temp)

            # Normalization
            # CPM, MAS5, fRMA, gcRMA, RMA, RMA+quantile normalization+FARMS, average, TMM,
            # RMA+quantile normalization, normexp

            ## Data is provided on a linear scale as Kallisto-estimated counts, i.e.  raw counts,
            ## convert it to TPM
            #df_data_temp = convertRawToTPM(df_data_temp, se_gene_length['hugo'])

            # Create columns as pandas.MultiIndex
            df_data_temp.columns = pd.MultiIndex.from_tuples([(dataset, sample) for sample in df_data_temp.columns.values])

            # Concatenate 'df_data_temp' to 'df_data'
            df_data = pd.concat([df_data, df_data_temp], sort=False, axis=1)

        if self.recordCollectionToBinaries:
            write(df_data, 'dev/df_data_R%s_%s' % (self.round, self.model_level))

        return df_data

    def run(self, num2keep=None):

        self.make_predictions(self.df_signature,
                            self.df_data,
                            self.model_name,
                            model_level=self.model_level,
                            df_gold_standard=self.df_gold_standard,
                            saveDir=self.outputDir,
                            round=self.round,
                            num2keep=num2keep)

        return

    def make_predictions(self, df_signature, df_data, modelName, model_level=None, df_gold_standard=None, saveDir=None, round=0, num2keep=None):

        df_signature.columns = df_signature.columns.get_level_values(-1)

        # Header of the predictions.csv
        predictions = [['dataset.name', 'sample.id', 'cell.type', 'prediction']]

        df_all_result = pd.DataFrame()

        for dataset in np.unique(df_data.columns.get_level_values(0).values):
            print('Evaluating dataset:', dataset)

            # Select current dataset and keep only genes that are present in signature matrix
            df_dataset = df_data.xs(dataset, level=0, axis=1).reindex(df_signature.index).fillna(0.)

            if not num2keep is None:
                print("selecting genes with a high coefficient of variance")
                top_coef_genes = self.keep_top_coef_of_var(df_dataset, df_signature, num2keep).unique()
                df_dataset = df_dataset.loc[top_coef_genes, :]
                df_signature = df_signature.loc[top_coef_genes, :]

            # Calculate fractions
            df_beta = self.get_beta_of_dataset(df_dataset, df_signature, modelName)

            for sample in df_beta.index:
                for celltype in df_beta.columns:
                    predictions.append([dataset, sample, celltype, df_beta.loc[sample, celltype]])

            # If the 'gold standard' is available make correlation plots calculate metrics.
            if not df_gold_standard is None:
                # Retrieve the gold standard for this dataset
                df_dataset_gold_standard = df_gold_standard.xs(key=dataset, level='dataset.name').unstack()

                # Reorder this dataset gold standard same as data and signature
                df_dataset_gold_standard = df_dataset_gold_standard.loc[df_dataset.columns][df_beta.columns]

                # Get plots and metrics for this dataset
                df_result = calculate_correlation_make_plots(df_beta.values, df_dataset_gold_standard.values, df_beta.columns, saveDir + str(dataset), noPlot=False)

                # Save metrics for this dataset
                df_temp_result = df_result.loc['Pearson'].to_frame().T
                df_temp_result.index = pd.Index([dataset])
                df_all_result = pd.concat([df_all_result, df_temp_result], axis=0, sort=False)

        if not df_gold_standard is None:
            print(self.model_name)
            print(df_all_result.T, '\n\n\n')
            df_all_result = df_all_result.unstack().dropna().to_frame().T
            values = df_all_result.unstack().dropna().to_frame().T.values[0]
            df_all_result.index = [np.round(np.mean(values[~np.isnan(values)]), 3)]
            df_all_result.to_excel('coef_of_var_params/df_all_result_%s_%s_R%s_marker_per_celltype_%s.xlsx'%(modelName, model_level, round, num2keep if not num2keep is None else ""))
            df_all_result.to_excel('df_all_result_%s_%s_R%s.xlsx'%(modelName, model_level, round))

        if not os.path.exists(os.path.join('output', '')):
            os.makedirs(os.path.join('output', ''))

        np.savetxt(os.path.join('output', 'predictions.csv'), np.vstack(predictions), fmt='%s', delimiter=',')

        return

    def keep_top_coef_of_var(self, df_data, df_signature, num2keep=30):
        # calculate the coefficient of variation of each gene in the dataset

        hard_coded = {"CD8.T.cells": 7} if self.model_level == "coarse" else {"naive.B.cells": 7}

        coef_of_var = (df_data.std(axis=1) / df_data.mean(axis=1)).abs()
        coef_of_var.sort_values(ascending=False, inplace=True) # sorted coefficient of variation

        keep_these = pd.Index([])
        for i, celltype in enumerate(df_signature.columns):
            genes = df_signature.loc[:, celltype]
            genes = (genes[genes == 1.]).index # markers for this cell type

            assert(not genes.duplicated().any())

            coef_for_celltype = coef_of_var.loc[genes] # CoV for this cell types markers
            coef_for_celltype.sort_values(ascending=False, inplace=True) # sort 'em

            if celltype in hard_coded.keys(): # fixing num for some celltypes
                high_var_genes = coef_for_celltype.iloc[:hard_coded[celltype]].index  # take the genes for this celltype with the highest CoV

            else:
                # high_var_genes = coef_for_celltype.iloc[:].index  # keep every gene
                high_var_genes = coef_for_celltype.iloc[:num2keep].index  # take the genes for this celltype with the highest CoV

            keep_these = keep_these.append(high_var_genes) # add those genes to the overall list...

        return keep_these

    def get_beta_of_dataset(self, df_data, df_signature, modelName, normalizationFunction = np.linalg.norm):

        '''Use signature matrix to deconvolve fractions (beta).
        df_data and df_signature should have identical index and no NaN values.
        Align them before calling this function

        normalizationFunction: function or None
            np.linalg.norm: divide every gene by its vector norm
            np.std: divide every gene by its standard deviation
            None: no normalization
        '''

        if not (df_data.index == df_signature.index).all():
            print('df_data and df_signature index is not aligned')

            return

        data = df_data.values.astype(float)
        signature = df_signature.values.astype(float)

        # Making sure that signature matrix has no weights, and no negatives
        signature[signature > 0.] = 1.
        signature[signature < 0.] = 0.

        # Normalize data gene-wise
        if not normalizationFunction is None:
            norm = normalizationFunction(data, axis=1)[:,None]
            data /= norm + 1. * (norm == 0.)

        # Score = (X1 + X2 + ...)/N
        if modelName == 'arithmetic_mean':

            beta = np.dot(data.T, signature) / data.shape[0]

        # Median of random samples
        elif modelName == 'median_random_samples':
            pick = lambda m, pool: np.random.choice(pool, size=min(m, len(pool)), replace=False)

            # The two approachaes below are statistically equivalent
            if False:
                # Sampling celltype by celltype, with smaller number of markers
                Q, M = 5000, 7
                beta = []
                for c in range(signature.shape[1]):
                    allBeta = data[np.vstack([pick(M, np.where(signature[:,c]==1)[0]) for i in range(Q)])].sum(axis=1)
                    allBeta = np.nanmedian(allBeta, axis=0)
                    print('*', end=' ', flush=True)
                    beta.append(allBeta)
                beta = np.vstack(beta).T
            else:
                # Sampling celltypes alltogether, with larger number of markers
                Q, M = 5000, 75
                idx = np.vstack([pick(M, list(range(len(signature)))) for i in range(Q)])
                allBeta = np.dstack([np.dot(data.T[:,idx[i,:]], signature[idx[i,:],:]) for i in range(Q)])
                allBeta[allBeta==0.] = np.nan
                beta = np.nanmedian(allBeta, axis=2)

        # Linear regression
        elif modelName == 'linear_regression':

            beta = linear_model.LinearRegression(fit_intercept=True, normalize=True).fit(signature, data).coef_

        # Elastic net regression
        elif modelName == 'elastic_net':

            beta = linear_model.ElasticNet(alpha=0.000001, l1_ratio=0.01, normalize=False).fit(signature, data).coef_

        # Arithmetic mean with penalty
        elif modelName == 'arithmetic_mean_minus':
            Z = scipy.stats.zscore(data.T, axis=1)

            beta = np.zeros((data.shape[1], signature.shape[1]))
            for j in range(signature.shape[1]):
                factor = np.zeros(data.T.shape) / data.shape[0]
                factor[((signature.T[j]==1.) * (Z<-0.1))] = 1.
                factor[((signature.T[j]==0.) * (Z>0.1))] = 1.
                factor = factor.sum(axis=1)

                beta[:,j] = np.dot(data.T, signature.T[j]) / (0.1 + factor)

        # Score = 2^((log2(X1) + log2(X2) + ...)/N)
        elif modelName == 'geometric_mean':
            beta = np.zeros((data.shape[1], signature.shape[1]))
            for i in range(data.shape[1]):
                for j in range(signature.shape[1]):
                    beta[i,j] = scipy.stats.mstats.gmean((data.T[i]*signature[:,j])[(data.T[i]*signature[:,j])>0.])

        # Score = (log2(X1) + log2(X2) + ...)/N, same as MCPCounter
        elif modelName == 'log_2_geometric_mean':
            beta = np.zeros((data.shape[1], signature.shape[1]))
            for i in range(data.shape[1]):
                for j in range(signature.shape[1]):
                    beta[i,j] = np.log2(scipy.stats.mstats.gmean((data.T[i]*signature[:,j])[(data.T[i]*signature[:,j])>0.]))

        else:
            print('Model type unknown')

            raise ValueError

        return pd.DataFrame(index=df_data.columns, columns=df_signature.columns, data=beta)
