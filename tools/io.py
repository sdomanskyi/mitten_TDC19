import os
import gzip
import pickle
import urllib.request
import tarfile
import pandas as pd
import numpy as np

def write(data, fileName):
    
    '''Pickle object into a file.'''

    try:
        with gzip.open(fileName + '.pklz','wb') as tempFile:
            pickle.dump(data, tempFile, protocol=4)

    except Exception as exception:
        print(exception)

    return

def read(fileName):

    '''Unpickle object from a file.'''

    if os.path.isfile(fileName + '.pklz'):
        try:
            with gzip.open(fileName + '.pklz','rb') as tempFile:
                data = pickle.load(tempFile)

                return data

        except Exception as exception:
            print(exception)

    return

def KeyInStore(key, file):

    try:
        with pd.HDFStore(file) as hdf5file:
            if "/" + key.strip("/") in hdf5file.keys():
                return True
            else:
                return False
    except Exception as exception:
        print(exception)

    return

def downloadFile(url, saveFolder):

    if not os.path.exists(saveFolder): 
        os.makedirs(saveFolder)
    
    fileName = url.strip('"').split('/')[-1:][0]
    path = os.path.join(saveFolder, fileName)

    if os.path.isfile(path):
        print('File has been downloaded already')
    else:
        print('Downloading file:', url.strip('"'), end='\t', flush=True)

        try:
            urllib.request.urlretrieve(url.strip('"'), path)
            print('Done', flush=True)

        except Exception as exception:
            print(exception)

    return

def loadData(name, model_level):

    if name == 'R1':
        fileName = 'dev/df_data_%s_%s' % (name, model_level)
        ds = read(fileName)

        if ds is None:
            # Read rounds data and record them to binaries, for quick access
            TumorDeconvolutionModels(model_level=model_level, input_name='input_r1.csv', round=1, recordCollectionToBinaries=True)
            ds = read(fileName)

        gs = pd.read_csv('input/gold_standards/lb_%s_r1.csv' % (model_level), index_col=[0,1,2], header=0)['measured']

    elif name == 'R2':
        fileName = 'dev/df_data_%s_%s' % (name, model_level)
        ds = read(fileName)
        
        if ds is None:
            # Read rounds data and record them to binaries, for quick access
            TumorDeconvolutionModels(model_level=model_level, input_name='input_r2.csv', round=2, recordCollectionToBinaries=True)
            ds = read(fileName)

        gs = pd.read_csv('input/gold_standards/lb_%s_r2.csv' % (model_level), index_col=[0,1,2], header=0)['measured']

    elif name == 'R3':
        fileName = 'dev/df_data_%s_%s' % (name, model_level)
        ds = read(fileName)
        
        if ds is None:
            # Read rounds data and record them to binaries, for quick access
            TumorDeconvolutionModels(model_level=model_level, input_name='input_r3_%s.csv' % (model_level), round=3, recordCollectionToBinaries=True)
            ds = read(fileName)

        gs = pd.read_csv('input/gold_standards/lb_%s_r3.csv' % (model_level), index_col=[0,1,2], header=0)['measured']

    elif name == 'DRC':
        ds = read('dev/df_samples_DRC_10by50_%s' % (model_level))
        gs = read('dev/df_gs_DRC_10by50_%s' % (model_level))['measured']

    elif name == 'SC':
        ds = read('dev/df_samples_SC_25by50_%s' % (model_level))
        gs = read('dev/df_gs_SC_25by50_%s' % (model_level))['measured']

    return ds, gs
