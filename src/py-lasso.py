# module
import argparse
import numpy as np
import pandas as pd
import os

from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import Lasso
from sklearn.multiclass import OneVsRestClassifier
from imblearn.over_sampling import SMOTE

# function
def feature_selection_LASSO(X, y, over_sampling=None):
    num_class = set(y)
    
    # SMOTE oversampling for minority
    if over_sampling:
      sm = SMOTE("minority",random_state=331)
      X, y = sm.fit_resample(X,y)
    
    if len(num_class) > 2:
        y = label_binarize(y, classes=list(num_class))

        # ML pipeline
        pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model', OneVsRestClassifier(Lasso(max_iter=99999)))])

        # grid search
        search = GridSearchCV(pipeline,
                      {'model__estimator__alpha':np.arange(0.01,3,0.001)},
                      cv = 10, scoring="neg_mean_squared_error",verbose=1
                      )

        search.fit(X,y)

        print(search.best_params_)
        coefficients = pd.DataFrame(search.best_estimator_.named_steps.model.coef_)
        return coefficients.abs().sum().to_list()

    else :
        pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model',Lasso(max_iter=99999))])

        # grid search
        search = GridSearchCV(pipeline,
                          {'model__alpha':np.arange(0.01,3,0.001)},
                          cv = 10, scoring="neg_mean_squared_error",verbose=1
                          )

        search.fit(X,y.values.ravel())

        print(search.best_params_)
        coefficients = search.best_estimator_.named_steps['model'].coef_
        return coefficients

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='Feature selection - LASSO')
    parser.add_argument('-b', '--base', required=True, type=str, help='Base direcotry')
    parser.add_argument('-x', '--feature', required=True, type=str, help='Feature - TXT file PATH')
    parser.add_argument('-y', '--predictor', required=True, type=str, help='Predictor variable - TXT file PATH')
    parser.add_argument('-o', '--oversampling', required=True, type=bool, help='Over sampling')
    static_args = parser.parse_args()

    # file path
    BASE = static_args.base
    X_PATH = static_args.feature
    Y_PATH = static_args.predictor
    OVER_SAMPLING = static_args.oversampling

    # run
    X = pd.read_csv(X_PATH)
    y = pd.read_csv(Y_PATH)

    rs = feature_selection_LASSO(X=X, y=y, over_sampling=OVER_SAMPLING)
    pd.DataFrame(rs, columns=['coef']).to_csv(BASE + "/lasso_result.csv")