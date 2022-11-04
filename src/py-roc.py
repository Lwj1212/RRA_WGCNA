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
def roc_acu_calculator(DF, feature_name, log_save, over_sampling, module_name):
  X = DF.iloc[:, 1:]
  y = DF.iloc[:, 0]
  
  # SMOTE oversampling for minority
  if over_sampling:
    sm = SMOTE("minority", random_state=331)
    X, y = sm.fit_resample(X,y)
    
  # multi-class detection
  num_class = set(y)
  if len(num_class) > 2:
    y = label_binarize(y, classes=list(num_class))
  
  # Learn to predict each class against the other 
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = .3, random_state=331)
  classifier = OneVsRestClassifier(SVC(kernel='linear', probability=True, random_state=331)) 
  y_score = classifier.fit(X_train, y_train).decision_function(X_test)
  
  if len(num_class) > 2:    
      
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(len(num_class)): 
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i]) 
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area 
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel()) 
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(num_class)):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= len(num_class)
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.4f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.4f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'darkgreen', 'violet', 'peru', 'gold'])
    for i, color in zip(range(len(num_class)), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (AUC={1:0.4f})'
                 ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(module_name + "_" + feature_name + " (Multi-class)")
    plt.legend(loc="lower right")
    plt.savefig(log_save + "/" + module_name + "_" + feature_name + '_ROC_AUC.png')
    plt.clf()
    
    # calculation
    y_prob = classifier.predict_proba(X_test)
    macro_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
    micro_roc_auc_ovr = roc_auc_score(
        y_test, y_prob, multi_class="ovr", average="micro")
    return {'macro' : macro_roc_auc_ovr, 'micro': micro_roc_auc_ovr}
    return macro_roc_auc_ovr
  else :
    fpr, tpr, thres = roc_curve(y_test, y_score)
    roc_auc = roc_auc_score(y_test, y_score)

    # roc curve
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.4f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(module_name + "_" + feature_name + " (Binary-class")
    plt.legend(loc="best")
    plt.savefig(log_save + "/" + module_name + "_" + feature_name + '_ROC_AUC.png')
    plt.clf()
    
    return roc_auc

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


    rs = roc_acu_calculator(DF=DF, feature_name=FEATURE_NAME, log_save=LOG_SAVE, over_sampling=OVER_SAMPLING, module_name=MODULE_NAME)
    
    # pd.DataFrame([rs], columns=['auc']).to_csv(BASE + "/ml_validation_result.csv")
    with open(BASE + '/ml_validation_result.json', 'w') as fp:
      json.dump(rs, fp)