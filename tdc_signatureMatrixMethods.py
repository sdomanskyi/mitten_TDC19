import os
from tools.io import *
import numpy as np
import pandas as pd
import scipy.stats
from tdc_models import TumorDeconvolutionModels

def getListOfExpressedGenes(names, model_level):

    fileName = os.path.join('dev','allGenes_%s' % ('_'.join(names)))

    if os.path.isfile(fileName + '.pklz'):
        try:
            return read(fileName)
        except Exception as exception:
            print(exception)

    allGenes = []

    for name in names:
        print('Collection:', name)

        ds = loadData(name, model_level)[0].fillna(0.)
            
        allGenes.extend(ds[np.abs(ds).sum(axis=1) > 0.].index.values.tolist())

    write(allGenes, fileName)

    return allGenes

def prepareSignatureMatrixZScore(df_training_data, groupLevel, nameSuffix, top = 50, z_cutoff = 0.67):

    #df_v = ((df_training_data.replace(0., np.nan).groupby(axis=1,level=groupLevel, sort=False).agg(np.nanstd).replace(np.nan, 0.) / \
    #    df_training_data.replace(0., np.nan).groupby(axis=1,level=groupLevel, sort=False).agg(np.nanmean).replace(np.nan, 0.)).replace(np.nan, 0.))

    #df_v = df_v.replace(0., 999.)
    #df_v += 1.

    # Group training data and get 'supporting' counts (counts of non-zero values divided by total
    # counts in the groups of samples)
    df_c = df_training_data.replace(0., np.nan).groupby(axis=1,level=groupLevel, sort=False).count().replace(np.nan, 0.) / \
        df_training_data.groupby(axis=1,level=groupLevel, sort=False).count()

    #  Group training data by median
    df_g = df_training_data.groupby(axis=1,level=groupLevel, sort=False).median()

    # Calculate z-scores
    df_z = pd.DataFrame(data=scipy.stats.zscore(df_g.values, axis=1), index=df_g.index, columns=df_g.columns)

    # Reduce those z-scores corresponding to lower 'supporting' counts
    df_z = df_c * df_z #/ df_v

    ## Drop genes expressed in multiple cell types
    #df_z = df_z.loc[[gene for gene in df_z.index if np.sort(df_z.loc[gene])[-2]<z_cutoff]]

    if type(df_z.columns) is pd.Index:
        df_z.columns = pd.MultiIndex.from_arrays([df_z.columns.values])
        df_g.columns = pd.MultiIndex.from_arrays([df_g.columns.values])

    # Drop cell types of 'Unknown', i.e.  cancer types
    try:
        df_z = df_z.drop(columns='unknown', level=0)
    except:
        pass

    # Drop genes that were markers in 'Unknown' columns
    df_z = df_z.loc[[gene for gene in df_z.index if max(df_z.loc[gene]) > z_cutoff]]

    tops = {column:top for column in df_z.columns}

    def selectGenes():

        # Keep only top-N genes with higest (z-score * count)
        selectedGenes = np.unique([item for sublist in [df_z.index[np.argsort(df_z[column].values)[-tops[column]:]].tolist() for column in df_z.columns] for item in sublist])

        temp = (df_z.loc[selectedGenes] > z_cutoff).sum(axis=0).values
        whereLow = np.argmin(temp)
        
        #print(whereLow, '\t', temp)

        return whereLow, selectedGenes

    df_z = df_z.loc[selectGenes()[1]]

    def func_to_0s_and_1s(subset):

        whereLarge = subset.values > z_cutoff
        subset[:] = 0.
        subset[whereLarge] = 1.

        return subset

    # Set values to 1 if the gene is a marker of that cell type and 0 otherwise
    df_z = df_z.apply(func_to_0s_and_1s, axis=1)

    #df_z *= df_g.loc[df_z.index][df_z.columns]

    return df_z

def filterSignatureMatrixCorr(model_level, L, backgroundGenesSet = None, signature_name = ''):
    
    def getGenesForCelltype(celltype, celltypes, model_level = 'fine', L = [], doTesting = False):

        def getBestCorrelatedGenes(names, backgroundGenesSet = None):

            df_signature_1 = pd.read_excel(signature_name, index_col=0, header=[0], skiprows=[1])
            backgroundGenesSet = df_signature_1[celltype][(df_signature_1[celltype] == 1.) & (df_signature_1.sum(axis=1) <= 3)].index.values

            df = pd.DataFrame()

            for name in names:
                ds, gs = loadData(name, model_level)

                for dataset in np.unique(ds.columns.get_level_values(0)):

                    ds_temp = ds.xs(level=0, key=dataset, axis=1).fillna(0.)
                    gs_temp = gs.xs(level=0, key=dataset, axis=0).unstack('sample.id')[ds_temp.columns]
                
                    try:
                        gs_temp = gs_temp.loc[celltype]
                    except Exception as exception:
                        print(exception)
                        continue

                    if np.isnan(gs_temp).all():
                        continue
                    else:
                        print('*', end=' ', flush=True) 
                        pass
                
                    if not backgroundGenesSet is None:
                        ds_temp = ds_temp.loc[ds_temp.index.intersection(np.unique(backgroundGenesSet))]

                    wh = np.where(~np.isnan(gs_temp.values))

                    c = np.array([np.corrcoef(ds_temp.values[i][wh], gs_temp.values[wh])[0,1] for i in range(len(ds_temp))])

                    df = pd.concat([df, pd.DataFrame(data=c[:,None], index=ds_temp.index, columns=[dataset])], axis=1, sort=False)

            se = df.min(axis=1, skipna=False).sort_values(ascending=False).dropna()

            if len(se) == 0.:
                print('Empty')

                return []

            se = se[se > (max(se) - 0.3)].iloc[:250]

            print(len(backgroundGenesSet), len(se), max(se), min(se))
            print(set(se.index.values))
            print()

            return set(se.index.values)
    
        def testCollections(selectedGenes, collections, model_level, celltype):

            def testDataset(name, selectedGenes):

                print('\n', name, end=':\t', flush=True)

                ds, gs = loadData(name, model_level)

                for dataset in np.unique(ds.columns.get_level_values(0)):

                    ds_temp = ds.xs(level=0, key=dataset, axis=1).fillna(0.)
                    gs_temp = gs.xs(level=0, key=dataset, axis=0).unstack('sample.id')[ds_temp.columns]

                    try:
                        gs_temp = gs_temp.loc[celltypes]
                    except Exception as exception:
                        print(exception)
                        continue

                    wh = np.where(~np.isnan(gs_temp.loc[celltype].values))[0]

                    if len(wh) > 0:
                        temp1 = ds_temp.loc[selectedGenes]
                        temp1 /= np.linalg.norm(temp1, axis=1)[:, None]
                    
                        temp1 = temp1.sum(axis=0)[wh]
                        temp2 = gs_temp.loc[celltype][wh]

                        corr = np.corrcoef(temp1, temp2)[0,1]

                        print('%s > %s' % (dataset, np.round(corr, 2)), end='\t', flush=True)

                return

            for collection in collections:
                testDataset(collection, selectedGenes)

            return

        print('\nCelltype:', celltype)

        selectedGenes = getBestCorrelatedGenes(L, backgroundGenesSet=backgroundGenesSet)
    
        if doTesting:
            testCollections(selectedGenes, ['R1', 'R2', 'R3'], model_level=model_level, celltype=celltype) # 'DRC', 'SC'

        return list(selectedGenes)

    if model_level == 'coarse':
        celltypes = ['B.cells', 
                        'CD4.T.cells', 
                        'CD8.T.cells', 
                        'NK.cells', 
                        'neutrophils', 
                        'monocytic.lineage', 
                        'fibroblasts', 
                        'endothelial.cells']
    elif model_level == 'fine':
        celltypes = ['NK.cells', 
                        'endothelial.cells', 
                        'fibroblasts', 
                        'macrophages', 
                        'memory.B.cells', 
                        'memory.CD4.T.cells', 
                        'memory.CD8.T.cells', 
                        'monocytes', 
                        'myeloid.dendritic.cells', 
                        'naive.B.cells', 
                        'naive.CD4.T.cells', 
                        'naive.CD8.T.cells', 
                        'neutrophils', 
                        'regulatory.T.cells']
    
    genesResultDict = {celltype: getGenesForCelltype(celltype, celltypes, model_level, L=L) for celltype in celltypes[:]}

    print('Combining markers into a DataFrame')
    df = pd.DataFrame(index=np.unique(np.hstack(list(genesResultDict.values()))), columns=celltypes).fillna(0.)

    try:
        for celltype in celltypes:
            df.loc[genesResultDict[celltype], celltype] = 1.
    except Exception as exception:
        print(exception)

    print('Saving')
    forIndex = [df.columns] if model_level == 'coarse' else [df.columns, df.columns]
    df.columns = pd.MultiIndex.from_arrays(forIndex)

    return df
