import os
import numpy as np
import pandas as pd
from tools.io import *
from tdc_signatureMatrixMethods import prepareSignatureMatrixZScore, filterSignatureMatrixCorr
from tdc_generateSyntheticData import generateSyntheticDatasets

celltypes = {'coarse': ['NK.cells', 'endothelial.cells', 'fibroblasts', 'B.cells', 'CD4.T.cells', 'CD8.T.cells', 'monocytic.lineage', 'neutrophils'],
             'fine': ['NK.cells','endothelial.cells','fibroblasts', 'macrophages', 'memory.B.cells', 'memory.CD4.T.cells', 'memory.CD8.T.cells', 'monocytes', 'myeloid.dendritic.cells', 'naive.B.cells', 'naive.CD4.T.cells', 'naive.CD8.T.cells', 'neutrophils', 'regulatory.T.cells']}

def loadBulkRNAseqData():

    dataFileName = 'data_v3.h5'

    # Coarse-grained names
    level_0_names = {'B cells' : 'B.cells',
                    'CD4 T cells' : 'CD4.T.cells',
                    'CD8 T cells' : 'CD8.T.cells',
                    'Endothelial cells' : 'endothelial.cells',
                    'Fibroblasts' : 'fibroblasts',
                    'Monocytic-lineage cells (i.e., monocytes, myeloid dendritic cells, macrophages)' : 'monocytic.lineage',
                    'Natural Killer cells' : 'NK.cells',
                    'Neutrophils' : 'neutrophils',
                    'Unknown' : 'unknown'}

    # Fine-grained names
    level_1_names = {'Breast cancer' : 'breast.cancer',
                    'Colorectal cancer' : 'colorectal.cancer',
                    'Endothelial cells' : 'endothelial.cells',
                    'Fibroblasts' : 'fibroblasts',
                    'Macrophages' : 'macrophages',
                    'Memory B cells' : 'memory.B.cells',
                    'Memory CD4 T cells' : 'memory.CD4.T.cells',
                    'Memory CD8 T cells' : 'memory.CD8.T.cells',
                    'Monocytes' : 'monocytes',
                    'Myeloid Dendritic cells' : 'myeloid.dendritic.cells',
                    'Naive B cells' : 'naive.B.cells',
                    'Naive CD4 T cells' : 'naive.CD4.T.cells',
                    'Naive CD8 T cells' : 'naive.CD8.T.cells',
                    'Natural Killer cells' : 'NK.cells',
                    'Neutrophils' : 'neutrophils',
                    'Regulatory T cells' : 'regulatory.T.cells'}

    if not KeyInStore('Subtypes_expression_by_samples', dataFileName):
        # Load data previously created by 'data_processing.py'
        df_bulkRNAseq = read('../Training Data/Subtypes_expression_by_samples')

        # Rename cell types and subtypes according to the DREAM Challenge names
        df_bulkRNAseq.columns = pd.MultiIndex.from_arrays([[level_0_names[item] for item in df_bulkRNAseq.columns.get_level_values(0).values],
                                    [level_1_names[item] for item in df_bulkRNAseq.columns.get_level_values(1).values],
                                    df_bulkRNAseq.columns.get_level_values(2).values,
                                    df_bulkRNAseq.columns.get_level_values(3).values])

        # Export data to HDF
        df_bulkRNAseq.to_hdf(dataFileName, key='Subtypes_expression_by_samples', mode='a', complevel=4, complib='zlib')
    else:
        print('Loading data from', dataFileName)
        df_bulkRNAseq = pd.read_hdf(dataFileName, key='Subtypes_expression_by_samples', mode='r')

    return df_bulkRNAseq

def loadScRNAseqData(model_level):

    if False:
        def coalesce(df):

            # Shuffle the cells
            df = df[np.random.choice(df.columns, size=len(df.columns), replace=False)]

            # Split cells into groups of 10
            dfs_split = np.array_split(df, np.int(df.shape[1] / 10.), axis=1)

            # Merge each group of cells of 10 into 1 by averaging
            dfs_grouped = [df_temp.groupby(axis=1, level=0).mean() for df_temp in dfs_split]

            # Combine pseudo-cells into one table
            dfs_merged = pd.concat(dfs_grouped, axis=1, sort=False, keys=range(len(dfs_grouped)))

            dfs_merged.columns = pd.MultiIndex.from_arrays([dfs_merged.columns.get_level_values(1), 
                                                            dfs_merged.columns.get_level_values(1), 
                                                            dfs_merged.columns.get_level_values(0)],
                                                            names=['celltype', 'celltype2', 'cell'])

            return dfs_merged

        if False:
            df_temp = pd.read_hdf('dev/coarseHCL.h5', key='key').fillna(0.)

            df_temp.xs(key='neutrophils', level=1, axis=1).to_hdf('dev/neutrophils.h5', key='key', mode='a', complevel=4, complib='zlib')
            df_temp.xs(key='fibroblasts', level=1, axis=1).to_hdf('dev/fibroblasts.h5', key='key', mode='a', complevel=4, complib='zlib')
            df_temp.xs(key='endothelial.cells', level=1, axis=1).to_hdf('dev/endothelial.cells.h5', key='key', mode='a', complevel=4, complib='zlib')

            df_temp = pd.concat([pd.read_hdf('dev/neutrophils.h5', key='key').fillna(0.),
                             pd.read_hdf('dev/fibroblasts.h5', key='key').fillna(0.),
                             pd.read_hdf('dev/endothelial.cells.h5', key='key').fillna(0.)], 
                            keys=['neutrophils', 'fibroblasts', 'endothelial.cells'],
                            axis=1, sort=False).fillna(0.)

            df_temp.columns = df_temp.columns.droplevel(1)
            df_temp = df_temp.loc[df_temp.sum(axis=1) > 0.]
            df_temp.columns.names = ['celltype', 'cell']

            df_temp.columns = pd.MultiIndex.from_arrays([df_temp.columns.get_level_values(0), 
                                                        df_temp.columns.get_level_values(0), 
                                                        df_temp.columns.get_level_values(1)],
                                                        names=['celltype', 'celltype2', 'cell'])

            df_temp.to_hdf('dev/neutro_fibro_endo.h5', key='key', mode='a', complevel=4, complib='zlib')

        else:
            df_neutro_fibro_endo = pd.read_hdf('dev/neutro_fibro_endo.h5', key='key')
            #27017 rows x 15000 columns
            df_neutro_fibro_endo = df_neutro_fibro_endo.apply(lambda q: q * 10000. / np.sum(q), axis=0)
            print(df_neutro_fibro_endo.sum(axis=0).head(10))

        if model_level == 'coarse':
            df = pd.read_hdf('dev/df_selected_%s.h5' % (model_level), key='df').fillna(0.)
            # 22754 rows x 16821 columns

            df.columns = pd.MultiIndex.from_arrays([df.columns.get_level_values(0), 
                                                    df.columns.get_level_values(0), 
                                                    df.columns.get_level_values(1)],
                                                    names=['celltype', 'celltype2', 'cell'])
   
        else:
            df = pd.read_hdf('dev/df_selected_%s.h5' % (model_level), key='df').fillna(0.)
            # 22754 rows x 22059 columns

            df.columns = pd.MultiIndex.from_arrays([df.columns.get_level_values(0), 
                                                    df.columns.get_level_values(0), 
                                                    df.columns.get_level_values(1)],
                                                    names=['celltype', 'celltype2', 'cell'])

        print('Converting to linear scale')
        df = np.power(2., df) - 1.
        df = df.apply(lambda q: q * 10000. / np.sum(q), axis=0)

        print(df.sum(axis=0).head(10))

        print('Cutting')
        for celltype in celltypes[model_level]:
            try:
                df1 = df_neutro_fibro_endo.xs(key=celltype, level=0, axis=1)
            except Exception as exception:
                print(exception)
                df1 = pd.DataFrame()
    
            try:
                df2 =  df.xs(key=celltype, level=0, axis=1)
            except Exception as exception:
                print(exception)
                df2 = pd.DataFrame()

            df_this_celltype = pd.concat([df1, df2], axis=1, sort=False).fillna(0.) 

            df_this_celltype = coalesce(df_this_celltype)
            write(df_this_celltype, 'dev/parts/%s/%s' % (model_level, celltype))
            print(df_this_celltype)

    df = pd.concat([read('dev/parts/%s/%s' % (model_level, celltype)) for celltype in celltypes[model_level]], axis=1, sort=False).fillna(0.)

    print(df)

    return df

np.random.seed(42)

if __name__ == '__main__':

    #The methods loadBulkRNAseqData and loadScRNAseqData must be overridden to use the latest data (the csv files, etc.)

    for model_level in ['coarse', 'fine']:
        for training_data_source in ['bulk', 'single cell']:
            if training_data_source == 'bulk':
                df_bulkRNAseq = loadBulkRNAseqData()

                generateSyntheticDatasets(df_bulkRNAseq, model_level, datasetNamePrefix='DRC', saveName='DRC', SC=False, numDatasets=10)

                print('Creating signature matrix Z-method (%s)' % (model_level))
                prepareSignatureMatrixZScore(df_bulkRNAseq, [0] if model_level == 'coarse' else [0, 1], model_level, top=100, z_cutoff=0.67).to_excel('Signature_matrix_%s_DRC.xlsx' % (model_level), index_label='index')

                print('Filtering signature matrix Corr-method (%s)' % (model_level))
                filterSignatureMatrixCorr(model_level, ['DRC'], signature_name='Signature_matrix_%s_DRC.xlsx' % (model_level)).to_excel('sig_corr_sel_%s_DRC_best.xlsx' % (model_level), index_label='index')

            if training_data_source == 'single cell':
                df_scRNAseq = loadScRNAseqData(model_level=model_level)

                generateSyntheticDatasets(df_scRNAseq, model_level, datasetNamePrefix='SC', saveName='SC', SC=True, numDatasets=25)

                print('Creating signature matrix Z-method (%s)' % (model_level))
                prepareSignatureMatrixZScore(df_scRNAseq, [0], model_level, top=50, z_cutoff=0.67).to_excel('Signature_matrix_%s_SC.xlsx' % (model_level), index_label='index')

                print('Filtering signature matrix Corr-method (%s)' % (model_level))
                filterSignatureMatrixCorr(model_level, ['SC'], signature_name='Signature_matrix_%s_SC.xlsx' % (model_level)).to_excel('sig_corr_sel_%s_SC_best.xlsx' % (model_level), index_label='index')

        # combo is the combination of bulk-derived and single-cell-derived signature matrices
        if True:
            celltypesSC = {'coarse': ['endothelial.cells', 'fibroblasts', 'CD4.T.cells', 'CD8.T.cells'],
                           'fine':   ['endothelial.cells', 'fibroblasts', 'memory.CD4.T.cells', 'memory.CD8.T.cells']}

            model_level = 'coarse'

            if model_level == 'coarse':
                df_DRC = pd.read_excel('sig_corr_sel_%s_DRC_best.xlsx' % (model_level), index_col=0, header=[0], skiprows=[1])
                df_SC = pd.read_excel('sig_corr_sel_%s_SC_best.xlsx' % (model_level), index_col=0, header=[0], skiprows=[1])
            else:
                df_DRC = pd.read_excel('sig_corr_sel_%s_DRC_best.xlsx' % (model_level), index_col=0, header=[0,1], skiprows=[2])
                df_SC = pd.read_excel('sig_corr_sel_%s_SC_best.xlsx' % (model_level), index_col=0, header=[0,1], skiprows=[2])

            df_DRC = df_DRC[df_DRC.columns[[(not item in celltypesSC[model_level]) for item in df_DRC.columns.get_level_values(0).values.tolist()]]]
            df_SC = df_SC[df_SC.columns[[(item in celltypesSC[model_level]) for item in df_SC.columns.get_level_values(0).values.tolist()]]]
            print(df_SC.shape, df_DRC.shape)

            df_DRC = df_DRC.loc[df_DRC.sum(axis=1)>0.]
            df_SC = df_SC.loc[df_SC.sum(axis=1)>0.]
            print(df_SC.shape, df_DRC.shape)

            df_combo = pd.concat([df_DRC, df_SC], axis=1).fillna(0.)
            print(df_combo.shape)

            df_combo = df_combo.loc[df_combo.sum(axis=1)>0.]
            print(df_combo.shape)

            df_combo = df_combo.sort_index()

            df_combo.to_excel('signature_matrix_combo_%s.xlsx' % (model_level), index_label='index')
