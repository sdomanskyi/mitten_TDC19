import os
import numpy as np
import pandas as pd
from tools.io import read, write

def generateSyntheticDatasets(df_data, groupLevel, numSamples = 50, numDatasets = 10, datasetNamePrefix = '', saveName = '', SC = True):

    try:
        df_data = df_data.drop(columns=['unknown'], level=0)
    except Exception as exception:
        print(exception)

    df_data = df_data.loc[(df_data.sum(axis=1) != 0).index]

    gLevel = 0 if groupLevel == 'coarse' else 1

    sizes = df_data.iloc[0:1].groupby(axis=1, level=gLevel, sort=False).count().T[df_data.index[0]].to_dict()
    print(sizes)

    gindex = df_data.iloc[0:1].groupby(axis=1, level=gLevel, sort=False).mean().columns

    print('Grouping index:', gindex)

    numMaxSamples = df_data.shape[1]

    df_samples = pd.DataFrame()
    df_gs = pd.DataFrame()

    for i in range(numDatasets):
        datasetName = datasetNamePrefix + '%s' % (i)

        # Randomly choose a subset of data samples
        selectedSamples = []
        for celltype in np.unique(df_data.columns.get_level_values(gLevel)):
            thisSamples = df_data.columns[df_data.columns.get_level_values(gLevel)==celltype]
            size = min(sizes[celltype], max(1 if SC else 1, int(np.random.rand()*sizes[celltype])))
            selectedSamples.extend(np.random.choice(thisSamples, size=size, replace=False))
        selectedSamples = pd.MultiIndex.from_tuples(selectedSamples)

        #Group the data by averaging the sapmples
        df_g = df_data[selectedSamples].groupby(axis=1, level=gLevel, sort=False).mean()

        # Create random fractions (a.k.a.  'gold standard')
        randomSequence = np.random.rand(numSamples, df_g.shape[1])**4.
        randomSequence /= np.sum(randomSequence, axis=1)[:, None]

        # Create random samples based on random fractions
        df_temp_random_samples = pd.DataFrame(data=np.dot(df_g.values, randomSequence.T), index=df_g.index, 
                                              columns=pd.MultiIndex.from_arrays([[datasetName] * numSamples, ['S%s' % (sample) for sample in range(numSamples)]]))

        df_samples = pd.concat([df_samples, df_temp_random_samples], sort=False, axis=1)

        temp_index = pd.MultiIndex.from_arrays([[datasetName] * numSamples * len(gindex), 
                                ['S%s' % (sample) for sample in range(numSamples) for j in range(len(gindex))], 
                                gindex.values.tolist() * numSamples], 
                                names=['dataset.name', 'sample.id', 'cell.type'])

        df_temp_gs = pd.DataFrame(index=pd.MultiIndex.from_arrays([[datasetName] * numSamples * len(df_g.columns), 
                                                                   ['S%s' % (sample) for sample in range(numSamples) for j in range(len(df_g.columns))], 
                                                                   df_g.columns.values.tolist() * numSamples], 
                                                                  names=['dataset.name', 'sample.id', 'cell.type']),
                                  columns=['measured'],
                                  data=np.hstack(list(randomSequence))).reindex(temp_index)

        df_gs = pd.concat([df_gs, df_temp_gs], sort=False, axis=0)

        print(datasetName, len(selectedSamples), len(df_g.columns), df_temp_random_samples.shape, df_temp_gs.shape)

    print(df_samples)
    print(df_gs)

    if not os.path.exists(os.path.join('dev','')):
        os.makedirs(os.path.join('dev',''))

    write(df_samples, os.path.join('dev', 'df_samples_%s_%sby%s_%s' % (saveName, numDatasets, numSamples, groupLevel)))
    write(df_gs, os.path.join('dev', 'df_gs_%s_%sby%s_%s' % (saveName, numDatasets, numSamples, groupLevel)))

    df_gs.to_csv(os.path.join('dev', 'df_gs_%s_%sby%s_%s.csv' % (saveName, numDatasets, numSamples, groupLevel)), na_rep='NA')

    return