#!/usr/bin/env python
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn import svm
from sklearn import metrics
import pylab as pl
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import pandas as pd
import numpy as np
from os.path import join
import itertools
from scipy import stats
import os
import math
import random
from imblearn.over_sampling import SMOTE
from collections import Counter


import pylab as pl
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

def plot_coefficients(coef, feature_names,path_affix, Features_showsize=10):

    positive_Features = min(len(feature_names[coef>0]),Features_showsize)
    negative_Features = min(len(feature_names[coef<0]),Features_showsize)
    top_positive_coefficients = np.argsort(coef)[-positive_Features:]
    top_negative_coefficients = np.argsort(coef)[:negative_Features]
    if (positive_Features>0 and negative_Features>0):
        top_coefficients = np.hstack([top_negative_coefficients, top_positive_coefficients])
    elif (positive_Features==0):
        top_coefficients = top_negative_coefficients
    else:
        top_coefficients = top_positive_coefficients
    # create plot
    plt.figure(figsize=(8, 5))
    colors = ['red' if c < 0 else 'blue' for c in coef[top_coefficients]]
    plt.bar(np.arange(positive_Features+negative_Features), coef[top_coefficients], color=colors)
    feature_names = np.array(feature_names)
    plt.xticks(np.arange(0.5, 0.5 + positive_Features+negative_Features), feature_names[top_coefficients], fontsize=14,rotation=60, ha='right')
    plt.savefig(path_affix+'_FeatureImportance.pdf',dpi=300)


def plot_roc(fpr, tpr,roc_auc, path_affix):

    pl.clf()
    pl.figure(figsize=(8, 5))
    pl.plot(fpr, tpr, label='ROC curve (area = %0.3f)' % roc_auc)
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('Receiver operating characteristic example')
    pl.legend(loc="lower right")
    pl.savefig(path_affix+'_roc.pdf',dpi=300)


def plot_prc(recall,precision,prc_auc, path_affix):

    pl.clf()
    pl.figure(figsize=(8, 5))
    pl.plot(recall, precision, label='PRC curve (area = %0.3f)' % prc_auc)
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.xlabel('Recall')
    pl.ylabel('Precision')
    pl.title('Precision recall characteristic example')
    pl.legend(loc="lower left")
    pl.savefig(path_affix+'_prc.pdf',dpi=300)


def plot_densComp(metric_name,metrics_random,metrics_biomarker, path_affix):

    df_metrics = pd.concat([pd.DataFrame({metric_name:(metrics_random[metric_name]).tolist(),
                                          'Features':np.repeat('random genes',metrics_random.shape[0]).tolist()}),
                            pd.DataFrame({metric_name:(metrics_biomarker[metric_name]).tolist(),
                                          'Features':np.repeat('biomarkers',metrics_biomarker.shape[0]).tolist()})])
    ks_statistic,ks_pvalue=stats.kstest(np.array(metrics_biomarker[metric_name]),
                                        np.array(metrics_random[metric_name]))
    ks_statistic,ks_pvalue=stats.kstest(np.array(metrics_biomarker[metric_name]),
                                        np.array(metrics_random[metric_name]))

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.text(0.8,1,'KS test (statistic = %0.2f, p-value = %0.2e)' %(ks_statistic,ks_pvalue))

    ax=sns.kdeplot(
       data=df_metrics, x=metric_name, hue="Features",
       fill=True, common_norm=False, palette="Paired",
       alpha=.5, linewidth=0,cut=0
    )
    plt.title('Density Plots of prediction performance across 100 repetitions')
    plt.xlabel('F1 score')
    plt.ylabel('Density')
    plt.savefig(path_affix+'_densityCompare.pdf',dpi=300)


def plot_boxComp(metric_name,metrics_random,metrics_biomarker, path_affix):

    df_metrics = pd.concat([pd.DataFrame({metric_name:(metrics_random[metric_name]).tolist(),
                                          'Features':np.repeat('random genes',metrics_random.shape[0]).tolist()}),
                            pd.DataFrame({metric_name:(metrics_biomarker[metric_name]).tolist(),
                                          'Features':np.repeat('biomarkers',metrics_biomarker.shape[0]).tolist()})])
    ks_statistic,ks_pvalue=stats.kstest(np.array(metrics_biomarker[metric_name]),
                                        np.array(metrics_random[metric_name]))

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.text(0.8,1,'KS test (statistic = %0.2f, p-value = %0.2e)' %(ks_statistic,ks_pvalue))

    ax=sns.kdeplot(
       data=df_metrics, x=metric_name, hue="Features",
       fill=True, common_norm=False, palette="Paired",
       alpha=.5, linewidth=0,cut=0
    )
    plt.title('Density Plots of prediction performance across 100 repetitions')
    plt.xlabel('F1 score')
    plt.ylabel('Density')
    plt.savefig(path_affix+'_densityCompare.pdf',dpi=300)


'''
Parameters
-----------
: a string join the path and filename affix for writing plots and CSVs
sampling_strategy: when dic,  pass to sample_strategy in SMOTE();
                   when float, pass to the ratio of data size in stratified bootstrapping
resampling: 'False', 'bootstrapping', 'smote',
            'auto'- if imbalanced classes do smote oversampling the minority class otherwise just stratified bootstraping with replacement to resample a same size of datast;

'''
def Ext_Val(X,y,feature_names,path_affix,sampling_strategy=1.0,resampling=False,repetitions=20,verbose=False):

    scaler = StandardScaler()
    scaler.fit(X)
    X_norm=scaler.transform(X)

    total_coefs =np.zeros(len(feature_names))
    total_probas_=[]
    total_y_test=[]
    total_y_pred = []
    all_accuracy=np.zeros(repetitions)
    all_recall=np.zeros(repetitions)
    all_precision=np.zeros(repetitions)
    all_fscore=np.zeros(repetitions)
    all_roc_auc=np.zeros(repetitions)
    all_prc_auc=np.zeros(repetitions)


    loo = LeaveOneOut()
    n_splits = loo.get_n_splits(X_norm)
    flag=0

    for iteration in range(0,repetitions):

        iteration_y_test=[]
        iteration_y_pred=[]
        iteration_probas_=[]
        iteration_flag=0
        for train_index, test_index in loo.split(X_norm):

            X_train, X_test = X_norm[train_index,:], X_norm[test_index,:]
            y_train, y_test = y[train_index], y[test_index]


            if resampling=='auto':
                # oversample the minority class to be equal of magority class by generating synthetic samples
                sm = SMOTE(sampling_strategy='all',k_neighbors=min(Counter(y_train)[0],Counter(y_train)[1])-1)
                #print('iteration'+str(iteration)+' orignial: ')
                #print(y_train)
                # stratified bootstrapping
                X_train, y_train = sm.fit_resample(X_train, y_train)

                #bootstrap_index = np.concatenate((np.random.choice(np.where(y_train==0)[0],size=Counter(y_train)[0]),np.random.choice(np.where(y_train==1)[0],size=Counter(y_train)[1])))
                #X_train = X_train[bootstrap_index,:]
                #y_train = y_train[bootstrap_index]
                #print('iteration'+str(iteration)+' after bootstrap: ')
                #print(y_train)

            elif resampling=='smote':
                sm = SMOTE(sampling_strategy=sampling_strategy,k_neighbors=min(Counter(y_train)[0],Counter(y_train)[1])-1)
                X_train, y_train = sm.fit_resample(X_train, y_train)
            elif resampling=='bootstrapping':
                bootstrap_index = np.concatenate((np.random.choice(np.where(y_train==0)[0],size=round(Counter(y_train)[0]*sampling_strategy)),np.random.choice(np.where(y_train==1)[0],size=round(Counter(y_train)[1]*sampling_strategy))))
                X_train = X_train[bootstrap_index,:]
                y_train = y_train[bootstrap_index]



            clf = svm.SVC(kernel='linear', probability=True,class_weight='balanced')
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            total_y_pred = np.concatenate((total_y_pred, y_pred))
            probas_ = clf.predict_proba(X_test)
            coef = clf.coef_.ravel()
            iteration_y_test = np.concatenate((iteration_y_test, y_test))
            total_y_test = np.concatenate((total_y_test, y_test))
            iteration_y_pred = np.concatenate((iteration_y_pred, y_pred))
            if flag==0:
                total_probas_=probas_
            else:
                total_probas_=np.concatenate((total_probas_, probas_))
                total_coefs += coef
            if iteration_flag==0:
                iteration_probas_=probas_
            else:
                iteration_probas_=np.concatenate((iteration_probas_, probas_))

            iteration_flag +=1
            flag+=1

            if verbose:
                print("True label:",y_test,"Predict label:",y_pred," Probobility:",probas_)

        accuracy = metrics.accuracy_score(iteration_y_test, iteration_y_pred)
        recall = metrics.recall_score(iteration_y_test, iteration_y_pred)
        precision = metrics.precision_score(iteration_y_test, iteration_y_pred)
        fscore = metrics.f1_score(iteration_y_test, iteration_y_pred)

        fprs, tprs, thresholds = metrics.roc_curve(iteration_y_test, iteration_probas_[:, 1])
        roc_auc = metrics.auc(fprs, tprs)
        precisions, recalls, thresholds = metrics.precision_recall_curve(iteration_y_test, iteration_probas_[:, 1])
        prc_auc = metrics.auc(recalls, precisions)

        all_accuracy[iteration] = accuracy
        all_recall[iteration] = recall
        all_precision[iteration] = precision
        all_fscore[iteration] = fscore
        all_prc_auc[iteration] = prc_auc
        all_roc_auc[iteration] = roc_auc

    #plot_coefficients(total_coefs/repetitions/n_splits, feature_names,path_affix=path_affix)

    print("confusion matrix:")
    conf_m = metrics.confusion_matrix(total_y_test, total_y_pred)/repetitions
    print(conf_m)

    avg_accuracy = np.mean(accuracy)
    avg_recall = np.mean(recall)
    avg_precision = np.mean(precision)
    avg_fscore = np.mean(fscore)


    fprs, tprs, thresholds = metrics.roc_curve(total_y_test, total_probas_[:, 1])
    roc_auc = metrics.auc(fprs, tprs)
    #plot_roc(fprs, tprs,roc_auc, path_affix)


    precisions, recalls, thresholds = metrics.precision_recall_curve(total_y_test, total_probas_[:, 1])
    prc_auc = metrics.auc(recalls, precisions)
    #plot_prc(recalls,precisions,prc_auc, path_affix)

    if verbose:
        print("fpr:", fpr, "\ntpr:",  tpr,"\nthresholds:", thresholds)
        print("precision:", precision, "\nrecall:", recall,"\nthresholds:", thresholds)


    df_metrics = pd.DataFrame({'accuracy': all_accuracy,'recall': all_recall,'precision': all_precision,"fscore":all_fscore,"roc_auc":all_roc_auc,"prc_auc":all_prc_auc})
    df_metrics.to_csv(path_affix+"_metrics.csv",index=True)

    df_metrics_statistic=pd.DataFrame({'std':df_metrics.std().tolist(),'max':df_metrics.max().tolist(),'mean':df_metrics.mean().tolist(),'min':df_metrics.min().tolist()})
    df_metrics_statistic.index=df_metrics.columns
    df_metrics_statistic.to_csv(path_affix+"_statistic.csv",index=True)

def Ext_Val_random(df_X,y,feature_size,path_affix,sampling_strategy=1.0,resampling=False,repetitions=20,verbose=False):



    all_accuracy=np.zeros(repetitions)
    all_recall=np.zeros(repetitions)
    all_precision=np.zeros(repetitions)
    all_fscore=np.zeros(repetitions)
    all_roc_auc=np.zeros(repetitions)
    all_prc_auc=np.zeros(repetitions)
    total_y_pred=[]
    total_y_test=[]

    for test_repetition in range(0,repetitions):

        locusID_random = random.sample(list(df_X.columns), feature_size)
        X = np.array(df_X[locusID_random])
        scaler = StandardScaler()
        scaler.fit(X)
        X_norm=scaler.transform(X)
        loo = LeaveOneOut()
        n_splits = loo.get_n_splits(X_norm)
        loo = LeaveOneOut()
        n_splits = loo.get_n_splits(X_norm)

        iteration_y_test=[]
        iteration_y_pred=[]
        iteration_probas_=[]

        flag=0

        for train_index, test_index in loo.split(X_norm):

            X_train, X_test = X_norm[train_index,:], X_norm[test_index,:]
            y_train, y_test = y[train_index], y[test_index]

            if resampling=='auto':
                # oversample the minority class to be equal of magority class by generating synthetic samples
                sm = SMOTE(sampling_strategy='all',k_neighbors=min(Counter(y_train)[0],Counter(y_train)[1])-1)
                # stratified bootstrapping
                X_train, y_train = sm.fit_resample(X_train, y_train)
                #bootstrap_index = np.concatenate((np.random.choice(np.where(y_train==0)[0],size=Counter(y_train)[0]),np.random.choice(np.where(y_train==1)[0],size=Counter(y_train)[1])))
                #X_train = X_train[bootstrap_index,:]
                #y_train = y_train[bootstrap_index]
            elif resampling=='smote':
                sm = SMOTE(sampling_strategy=sampling_strategy,k_neighbors=min(Counter(y_train)[0],Counter(y_train)[1])-1)
                X_train, y_train = sm.fit_resample(X_train, y_train)
            elif resampling=='bootstrapping':
                bootstrap_index = np.concatenate((np.random.choice(np.where(y_train==0)[0],size=round(Counter(y_train)[0]*sampling_strategy)),np.random.choice(np.where(y_train==1)[0],size=round(Counter(y_train)[1]*sampling_strategy))))
                X_train = X_train[bootstrap_index,:]
                y_train = y_train[bootstrap_index]


            clf = svm.SVC(kernel='linear', probability=True,class_weight='balanced')
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            probas_ = clf.predict_proba(X_test)
            iteration_y_test=np.concatenate((iteration_y_test, y_test))
            iteration_y_pred=np.concatenate((iteration_y_pred, y_pred))
            total_y_test = np.concatenate((total_y_test, y_test))
            total_y_pred = np.concatenate((total_y_pred, y_pred))

            if flag==0:
                iteration_probas_=probas_
            else:
                iteration_probas_=np.concatenate((iteration_probas_, probas_))

            flag +=1
            if verbose:
                print("True label:",y_test,"Predict label:",y_pred," Probobility:",probas_)




        accuracy = metrics.accuracy_score(iteration_y_test, iteration_y_pred)
        recall = metrics.recall_score(iteration_y_test, iteration_y_pred)
        precision = metrics.precision_score(iteration_y_test, iteration_y_pred)
        fscore = metrics.f1_score(iteration_y_test, iteration_y_pred)

        fprs, tprs, thresholds = metrics.roc_curve(iteration_y_test, iteration_probas_[:, 1])
        roc_auc = metrics.auc(fprs, tprs)
       # plot_roc(fpr, tpr,roc_auc, path_affix)
        precisions, recalls, thresholds = metrics.precision_recall_curve(iteration_y_test, iteration_probas_[:, 1])
        prc_auc = metrics.auc(recalls, precisions)
        #plot_prc(recall,precision,prc_auc, path_affix)

        all_accuracy[test_repetition] = accuracy
        all_recall[test_repetition] = recall
        all_precision[test_repetition] = precision
        all_fscore[test_repetition] = fscore
        all_prc_auc[test_repetition] = prc_auc
        all_roc_auc[test_repetition] = roc_auc

        if verbose:
            print("fpr:", fpr, "\ntpr:",  tpr,"\nthresholds:", thresholds)
            print("precision:", precision, "\nrecall:", recall,"\nthresholds:", thresholds)

    print("confusion matrix:")
    conf_m = metrics.confusion_matrix(total_y_test, total_y_pred)/repetitions
    print(conf_m)
    df_metrics = pd.DataFrame({'accuracy': all_accuracy,'recall': all_recall,'precision': all_precision,"fscore":all_fscore,"roc_auc":all_roc_auc,"prc_auc":all_prc_auc})
    df_metrics.to_csv(path_affix+"_metrics.csv",index=True)


    df_metrics_statistic=pd.DataFrame({'std':df_metrics.std().tolist(),'max':df_metrics.max().tolist(),'mean':df_metrics.mean().tolist(),'min':df_metrics.min().tolist()})
    df_metrics_statistic.index=df_metrics.columns
    df_metrics_statistic.to_csv(path_affix+"_statistic.csv",index=True)
