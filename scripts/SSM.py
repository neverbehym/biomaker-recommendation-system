import pandas as pd
import numpy as np
from os.path import join
import itertools
from scipy import stats
import os
import math
import random
import pylab as pl
import matplotlib.pyplot as plt
import seaborn as sns
from util_ssm import *
import warnings
warnings.filterwarnings('ignore')

Panellog_FILE = str(sys.argv[1])
Panel_ID= int(sys.argv[2])
### read the dataframe of panel log with the index being panel_ID and column 'genes' being biomarker gene names seperated with ";"
df_panels = pd.read_csv(Panellog_FILE,index_col=0)

Genelist_FILE = str(sys.argv[3])
### read the dataframe of gene list with the index being locus_tag and column 'name' being gene names
df_genes = pd.read_csv(Genelist_FILE,index_col=0)

panel_size = df_panels.loc[Panel_ID,'size']
panel_genes = df_panels.loc[Panel_ID,'genes'].split(sep=';')
panel_locus = df_genes.loc[[x in panel_genes for x in df_genes['name']]].index.to_numpy()

data_DIR = str(sys.argv[4])
filenames = os.listdir(data_DIR)

biomarkerORrandom = str(sys.argv[5])
if biomarkerORrandom == 'random':
    result_DIR = '../results/SSM/RandomPanels/panel'+str(Panel_ID)
    panel_locus =  random.sample(list(df_genes.index), panel_size)
    panel_genes = df_genes.loc[panel_locus,'name']
else:
    result_DIR = '../results/SSM/BiomarkerPanels/panel'+str(Panel_ID)

if(os.path.isdir(result_DIR)==False):
    os.mkdir(makedirs)

for filename in filenames:

    ### token indicating the data condition
    data_token = filename.split(sep='.')[0]
    ### full file path to read the data
    data_FILE = join(data_DIR,filename)
    ### read data contains gene expression profiles in treatment samples and control samples
    df_data = pd.read_csv(data_FILE,index_col=0).transpose()
    y = df_data['label']

    ### get reduced dataset
    panel_locus_cur = [x for x in panel_locus if x in df_data.columns]
    panel_genes_cur = np.array(df_genes.loc[panel_locus_cur,'name'])
    print(str(len(panel_locus_cur))+" biomarkers are measured in the"+ data_token+" dataset.")
    X = df_data[panel_locus_cur]

    print (data_token+' - panel'+str(Panel_ID))
    Ext_Val(X,y,panel_genes_cur,join(result_DIR,data_token),resampling='auto',repetitions=3)
