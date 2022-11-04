# module
import argparse
import numpy as np
from numpy import interp
import pandas as pd
import json
from itertools import cycle

from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE

# function
def roc_acu_calculator_cv(DF, feature_name, log_save, over_sampling, module_name):
    X = DF.iloc[:, 1:]
    y = DF.iloc[:, 0]

    # SMOTE oversampling for minority
    if over_sampling:
        sm = SMOTE(sampling_strategy="minority", random_state=331)
        X, y = sm.fit_resample(X,y)

    # multi-class detection
    num_class = set(y)
    if len(num_class) > 2:
        y = label_binarize(y, classes=list(num_class))

    fold_cnt = 5
    kfold = KFold(n_splits=fold_cnt)
    scores = []
    pipe_line = make_pipeline(StandardScaler(),
                             OneVsRestClassifier(SVC(kernel='linear', probability=True, random_state=331)))

    roc_auc_fold = dict()
    fpr_fold = dict()
    tpr_fold = dict()

    for k, (train, test) in enumerate(kfold.split(DF)):
        y_score = pipe_line.fit(X.iloc[train, :], y[train]).decision_function(X.iloc[test, :])

        if len(num_class) > 2:
            fpr, tpr, roc_auc = roc_auc_function(y_test=y[test], y_score=y_score, num_class=num_class)
            roc_auc_fold[k] = roc_auc['micro']
            fpr_fold[k] = fpr['micro']
            tpr_fold[k] = tpr['micro']

        else :
            fpr, tpr, thres = roc_curve(y[test], y_score, drop_intermediate=False)
            roc_auc_fold[k] = roc_auc_score(y[test], y_score)
            fpr_fold[k] = fpr
            tpr_fold[k] = tpr


    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr_fold[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    # for i in range(len(num_class)):
    for i in range(fold_cnt):
        mean_tpr += interp(all_fpr, fpr_fold[i], tpr_fold[i])

    # Finally average it and compute AUC
    # mean_tpr /= len(num_class)
    mean_tpr /= fold_cnt
    fpr_fold["macro"] = all_fpr
    tpr_fold["macro"] = mean_tpr
    roc_auc_fold["macro"] = auc(fpr_fold["macro"], tpr_fold["macro"])

    # Plot all ROC curves
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr_fold["macro"], tpr_fold["macro"],
             label='macro-average ROC curve (AUC = {0:0.3f})'
                   ''.format(roc_auc_fold["macro"]),
             color='red', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'darkgreen', 'violet', 'peru', 'gold'])
    for i, color in zip(range(fold_cnt), colors):
        plt.plot(fpr_fold[i], tpr_fold[i], color=color, lw=lw, alpha=0.2,
                 label='ROC curve fold-{0} (AUC={1:0.3f})'
                 ''.format(i + 1, roc_auc_fold[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(module_name + "_" + feature_name + "_Multi-class-CV")
    plt.legend(loc="lower right")
    plt.savefig(log_save + "/Step4_" + module_name + "_" + feature_name + '_ROC_AUC_CV.png')
    plt.close()

    return roc_auc_fold["macro"]

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='ROC curve')
    parser.add_argument('-b', '--base', required=True, type=str, help='Base direcotry')
    parser.add_argument('-d', '--dataframe', required=True, type=str, help='DataFrame')
    parser.add_argument('-l', '--logsave', required=True, type=str, help='ROC curve plot save path')
    parser.add_argument('-f', '--featurename', required=True, type=str, help='ROC curve plot - title')
    parser.add_argument('-m', '--moduleename', required=True, type=str, help='ROC curve plot - title')
    parser.add_argument('-o', '--oversampling', required=True, type=bool, help='Over sampling')
    static_args = parser.parse_args()

    # file path
    BASE = static_args.base
    DATA_FRAME = static_args.dataframe
    LOG_SAVE = static_args.logsave
    FEATURE_NAME = static_args.featurename
    MODULE_NAME = static_args.moduleename
    OVER_SAMPLING = static_args.oversampling

    # run
    DF = pd.read_csv(DATA_FRAME)

    rs = roc_acu_calculator_cv(DF=DF, feature_name=FEATURE_NAME, log_save=LOG_SAVE, over_sampling=OVER_SAMPLING, module_name=MODULE_NAME)
        
    # pd.DataFrame([rs], columns=['auc']).to_csv(BASE + "/ml_validation_result.csv")
    with open(BASE + '/ml_validation_result.json', 'w') as fp:
      json.dump(rs, fp)
    