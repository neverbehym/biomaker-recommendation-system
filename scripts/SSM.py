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
from util_val import *
import warnings
warnings.filterwarnings('ignore')


DATA_DIR = join('./','valdata')
PANEL_DIR = join('./','valdata/BiomarkerPanels949/panel')
RESULT_DIR =join('./','results/validation949')
gene_list = pd.read_csv(join(DATA_DIR,'gene_completelist.csv'),index_col=0)


# read data
potassium_log2intensity = pd.read_csv(join(DATA_DIR,'potassium_log2intensity.csv'),index_col=0).transpose()
y_potassium = potassium_log2intensity['potassium_labels']

Mascher_log2tpm = pd.read_csv(join(DATA_DIR,'Mascher_log2tpm.csv'),index_col=0).transpose()
y_antibiotic = Mascher_log2tpm['antibiotic_labels']

HU_log2rpkm = pd.read_csv(join(DATA_DIR,'HU_log2rpkm.csv'),index_col=0).transpose()
y_HU = HU_log2rpkm['HU_labels']

Rath_log2intensity = pd.read_csv(join(DATA_DIR,'Rath_log2intensity.csv'),index_col=0).transpose()
y_salinity = Rath_log2intensity['Salinity_labels']
y_GB = Rath_log2intensity['GB_labels']

Heat_log2rpkm = pd.read_csv(join(DATA_DIR,'Heat_log2rpkm.csv'),index_col=0).transpose()
y_Heat = Heat_log2rpkm['Heat_labels']

Cold_log2rpkm = pd.read_csv(join(DATA_DIR,'Cold_log2rpkm.csv'),index_col=0).transpose()
y_Cold = Cold_log2rpkm['Cold_labels']

H2O2_log2tpm = pd.read_csv(join(DATA_DIR,'H2O2_log2tpm.csv'),index_col=0).transpose()
y_H2O2 = H2O2_log2tpm['H2O2_labels']

pressure_log2tpm = pd.read_csv(join(DATA_DIR,'pressure_log2tpm.csv'),index_col=0).transpose()
y_pressure = pressure_log2tpm['Pressure_labels']

sta_relintensity = pd.read_csv(join(DATA_DIR,'sta_relintensity.csv'),index_col=0).transpose()
y_sta =sta_relintensity['Stationary_labels']

deeps_log2tpm = pd.read_csv(join(DATA_DIR,'deeps_log2tpm.csv'),index_col=0).transpose()
y_deeps = deeps_log2tpm['DeepStarvation_labels']

for i in range(5,11):    #<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<---<-<-<-<-<-<---<-<-<-<-<-<---<-<-<-<-<-<---
    if(os.path.isdir(RESULT_DIR+'/panelRandomSize'+str(i))==False):
        os.mkdir(RESULT_DIR+'/panelRandomSize'+str(i))
    print ('antibiotic - random panels')
    Ext_Val_random(Mascher_log2tpm,y_antibiotic,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/antibiotic',resampling='auto',repetitions=200,verbose=False)
    print ('potassium - random panels')
    Ext_Val_random(potassium_log2intensity,y_potassium,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/potassium',resampling='auto',repetitions=200,verbose=False)
    print ('HU - random panels')
    Ext_Val_random(HU_log2rpkm,y_HU,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/HU',resampling='auto',repetitions=200,verbose=False)
    print ('salinity - random panels')
    Ext_Val_random(Rath_log2intensity,y_salinity,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/salinity',resampling='auto',repetitions=200,verbose=False)
    print ('GB - random panels')
    Ext_Val_random(Rath_log2intensity,y_GB,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/GB',resampling='auto',repetitions=200,verbose=False)
    print ('heat - random panels')
    Ext_Val_random(Heat_log2rpkm,y_Heat,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/Heat',resampling='auto',repetitions=200,verbose=False)
    print ('cold - random panels')
    Ext_Val_random(Cold_log2rpkm,y_Cold,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/Cold',resampling='auto',repetitions=200,verbose=False)
    print ('H2O2 - random panels')
    Ext_Val_random(H2O2_log2tpm,y_H2O2,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/H2O2',resampling='auto',repetitions=200,verbose=False)
    print ('pressure - random panels')
    Ext_Val_random(pressure_log2tpm,y_pressure,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/pressure',resampling='auto',repetitions=200,verbose=False)
    print ('stationary - random panels')
    Ext_Val_random(sta_relintensity,y_sta,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/stationary',resampling='auto',repetitions=200,verbose=False)
    print ('deep starvation - random panels')
    Ext_Val_random(deeps_log2tpm,y_deeps,i,RESULT_DIR+'/panelRandomSize'+str(i)+'/DeepStarvation',resampling='auto',repetitions=200,verbose=False)

#
#
# for panelID in range(1,950): #<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<-<---<-<-<-<-<-<---<-<-<-<-<-<---<-<-<-<-<-<---
#     ##print('panel'+str(panelID))
#     with open(PANEL_DIR+str(panelID)+'.txt') as f:
#         panel = [line.rstrip() for line in f]
#     row_mask=[x in panel for x in gene_list['name']]
#     panel_locusID = gene_list.loc[row_mask].index.to_numpy()
#     gene_list.loc[panel_locusID]
#     RESULT_DIR_cur = RESULT_DIR+'/panel'+str(panelID)
#     if(os.path.isdir(RESULT_DIR_cur)==False):
#         os.mkdir(RESULT_DIR_cur)
#
#     #################################### antibiotic ############################
#     # get reduced dataset
#     Mascher_locusID_biomarkers = [x for x in panel_locusID if x in Mascher_log2tpm.columns]
#     #print(str(len(Mascher_locusID_biomarkers))+" biomarkers are measured in the dataset for antibiotic condition.")
#
#     Mascher_genename_reduced = np.array(gene_list.loc[Mascher_locusID_biomarkers,'name'])
#     Mascher_log2tpm_reduced=Mascher_log2tpm[Mascher_locusID_biomarkers]
#
#     print ('antibiotic - panel'+str(panelID))
#     Ext_Val(Mascher_log2tpm_reduced,y_antibiotic,Mascher_genename_reduced,RESULT_DIR_cur+'/antibiotic',resampling='auto',repetitions=3)
#     Mascher_metrics = pd.read_csv(RESULT_DIR_cur+'/antibiotic_metrics.csv',index_col=0)
#
# #################################### potassium ############################
#     # get reduced datasets
#     potassium_locusID_biomarkers = [x for x in panel_locusID if x in potassium_log2intensity.columns]
#     #print(str(len(potassium_locusID_biomarkers))+" biomarkers are measured in the dataset for salinity and Glycine betaine conditions")
#     potassium_genename_reduced = np.array(gene_list.loc[potassium_locusID_biomarkers,'name'])
#     potassium_log2intensity_reduced=potassium_log2intensity[potassium_locusID_biomarkers]
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('potassium - panel'+str(panelID))
#     Ext_Val(potassium_log2intensity_reduced,y_potassium,potassium_genename_reduced,RESULT_DIR_cur+'/potassium',resampling='auto',repetitions=3)
#
# #################################### HU ############################
#
#     # get reduced datasets
#     HU_locusID_biomarkers = [x for x in panel_locusID if x in HU_log2rpkm.columns]
#     #print(str(len(Rath_locusID_biomarkers))+" biomarkers are measured in the dataset for salinity and Glycine betaine conditions")
#     HU_genename_reduced = np.array(gene_list.loc[HU_locusID_biomarkers,'name'])
#     HU_log2rpkm_reduced=HU_log2rpkm[HU_locusID_biomarkers]
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('HU - panel'+str(panelID))
#     Ext_Val(HU_log2rpkm_reduced,y_HU,HU_genename_reduced,RESULT_DIR_cur+'/HU',resampling='auto',repetitions=3)
#
# ####################################salinity############################
#     # get reduced datasets
#     Rath_locusID_biomarkers = [x for x in panel_locusID if x in Rath_log2intensity.columns]
#     #print(str(len(Rath_locusID_biomarkers))+" biomarkers are measured in the dataset for salinity and Glycine betaine conditions")
#     Rath_genename_reduced = np.array(gene_list.loc[Rath_locusID_biomarkers,'name'])
#     Rath_log2intensity_reduced=Rath_log2intensity[Rath_locusID_biomarkers]
#
#     print ('salinity - panel'+str(panelID))
#     Ext_Val(Rath_log2intensity_reduced,y_salinity,Rath_genename_reduced,RESULT_DIR_cur+'/salinity',resampling='auto',repetitions=3)
#     salinity_metrics = pd.read_csv(RESULT_DIR_cur+'/salinity_metrics.csv',index_col=0)
#
#   ####################################Glycine betaine############################
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('GB - panel'+str(panelID))
#     Ext_Val(Rath_log2intensity_reduced,y_GB,Rath_genename_reduced,RESULT_DIR_cur+'/GB',resampling='auto',repetitions=3,verbose=False)
#     GB_metrics = pd.read_csv(RESULT_DIR_cur+'/GB_metrics.csv',index_col=0)
#
#   #################################### Heat stroke ############################
#   # get reduced dataset
#     Heat_locusID_biomarkers = [x for x in panel_locusID if x in Heat_log2rpkm.columns]
#     #print(str(len(Heat_locusID_biomarkers))+" biomarkers are measured in the dataset for heat condition.")
#
#     Heat_genename_reduced = np.array(gene_list.loc[Heat_locusID_biomarkers,'name'])
#     Heat_log2rpkm_reduced=Heat_log2rpkm[Heat_locusID_biomarkers]
#   ## cross-validation on biomarkers, 100 repetitions
#     print ('heat - panel'+str(panelID))
#     Ext_Val(Heat_log2rpkm_reduced,y_Heat,Heat_genename_reduced,RESULT_DIR_cur+'/Heat',resampling='auto',repetitions=3)
#     Heat_metrics = pd.read_csv(RESULT_DIR_cur+'/Heat_metrics.csv',index_col=0)
#
#
#
#
#     #################################### Cold shock ############################
#     # get reduced dataset
#     Cold_biomarkers = [x for x in panel_locusID if x in Cold_log2rpkm.columns]
#     #print(str(len(Cold_biomarkers))+" biomarkers are measured in the dataset for cold condition.")
#     Cold_genename_reduced = np.array(gene_list.loc[Cold_biomarkers,'name'])
#     Cold_log2rpkm_reduced=Cold_log2rpkm[Cold_biomarkers]
#
#
#     print ('cold - panel'+str(panelID))
#     Ext_Val(Cold_log2rpkm_reduced,y_Cold,Cold_genename_reduced,RESULT_DIR_cur+'/Cold',resampling='auto',repetitions=3)
#
#
#
#   #################################### Oxidative stress ############################
#     # get reduced dataset
#     H2O2_biomarkers = [x for x in panel_locusID if x in H2O2_log2tpm.columns]
#     #print(str(len(H2O2_biomarkers))+" biomarkers are measured in the dataset for H2O2 condition.")
#
#     H2O2_genename_reduced = np.array(gene_list.loc[H2O2_biomarkers,'name'])
#     H2O2_log2tpm_reduced=H2O2_log2tpm[H2O2_biomarkers]
#
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('H2O2 - panel'+str(panelID))
#     Ext_Val(H2O2_log2tpm_reduced,y_H2O2,H2O2_genename_reduced,RESULT_DIR_cur+'/H2O2',resampling='auto',repetitions=3)
#
#
#     #################################### Pressure stress ############################
#     # get reduced dataset
#     pressure_locusID_biomarkers = [x for x in panel_locusID if x in pressure_log2tpm.columns]
#     #print(str(len(pressure_locusID_biomarkers))+" biomarkers are measured in the dataset for pressure condition.")
#
#     pressure_genename_reduced = np.array(gene_list.loc[pressure_locusID_biomarkers,'name'])
#     pressure_log2tpm_reduced=pressure_log2tpm[pressure_locusID_biomarkers]
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('pressure - panel'+str(panelID))
#     Ext_Val(pressure_log2tpm_reduced,y_pressure,pressure_genename_reduced,RESULT_DIR_cur+'/pressure',resampling='auto',repetitions=3)
#
#
#     #################################### Staionary Phase ############################
#     # get reduced dataset
#     sta_locusID_biomarkers = [x for x in panel_locusID if x in sta_relintensity.columns]
#     #print(str(len(sta_locusID_biomarkers))+" biomarkers are measured in the dataset for stationary condition.")
#     sta_genename_reduced = np.array(gene_list.loc[sta_locusID_biomarkers,'name'])
#     sta_relintensity_reduced=sta_relintensity[sta_locusID_biomarkers]
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('stationary - panel'+str(panelID))
#
#     Ext_Val(sta_relintensity_reduced,y_sta,sta_genename_reduced,RESULT_DIR_cur+'/stationary',resampling='auto',repetitions=3)
#     sta_metrics = pd.read_csv(RESULT_DIR_cur+'/stationary_metrics.csv',index_col=0)
#
#     #################################### Deep Starvation############################
#     # get reduced dataset
#     deeps_locusID_biomarkers = [x for x in panel_locusID if x in deeps_log2tpm.columns]
#     #print(str(len(deeps_locusID_biomarkers))+" biomarkers are measured in the dataset for deep starvation condition.")
#
#     deeps_genename_reduced = np.array(gene_list.loc[deeps_locusID_biomarkers,'name'])
#     deeps_log2tpm_reduced=deeps_log2tpm[deeps_locusID_biomarkers]
#
#     ## cross-validation on biomarkers, 100 repetitions
#     print ('deep starvation  - panel'+str(panelID))
#     Ext_Val(deeps_log2tpm_reduced,y_deeps,deeps_genename_reduced,RESULT_DIR_cur+'/DeepStarvation',resampling='auto',repetitions=3)
#     deeps_metrics = pd.read_csv(RESULT_DIR_cur+'/DeepStarvation_metrics.csv',index_col=0)
