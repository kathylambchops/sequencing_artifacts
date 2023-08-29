#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import os.path
from collections import Counter
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.impute import SimpleImputer 
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer

import time
import joblib
from sklearn.model_selection import cross_val_score

import get_metrics
import lazypredict
from lazypredict.Supervised import LazyClassifier

import sklearn
import xgboost
import lightgbm

from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from random import shuffle

#df = pd.read_pickle('ML_df.pkl')
df = pd.read_pickle('ML_df_joined_LR_feats.pkl')
#df = pd.read_pickle('ML_df_add_LR_feats.pkl')

X = df.drop(["IS_ARTIFACT"], axis=1) # drop labels
y = df["IS_ARTIFACT"].copy() # labels

cat_attribs = ['REF','ALT'] 
X_train_num = X.drop(cat_attribs, axis=1) 
num_attribs = list(X_train_num)

#all_features = num_attribs + cat_attribs

def select_random_features(num):
    all_features = num_attribs + cat_attribs
    shuffle(all_features)
    return all_features[:num]

# Investigate predictive power of features
# Preprocesses dataframe with a certain set of features to keep
def preprocess_with_features(df, features_to_keep):
    num_pipeline =  Pipeline([
    ('imputer', SimpleImputer(strategy="median")), # fill in missing values with median values
    ('std_scaler', StandardScaler()),              # feature scaling
        ])
    
    num_features_to_use =   list(set(num_attribs)   & set(features_to_keep))
    cat_features_to_use = list(set(cat_attribs) & set(features_to_keep))
    
    # Randomly selecting features. If by chance, we don't select any numerical (or categorical) features, 
    # we shouldn't should use the numerical part of the pipeline
    l = []
    if len(num_features_to_use) > 0:
        l.append(("num", num_pipeline, num_features_to_use))
    if len(cat_features_to_use) > 0:
        l.append(("cat", OneHotEncoder(), cat_features_to_use))
    full_pipeline = ColumnTransformer(l)
    

    # df_transformed is the full feature set. Now ready to use.
    #fitted = full_pipeline.fit(df)
    #df_transformed = fitted.transform(df)
    
    X = df.drop(["IS_ARTIFACT"], axis=1) # drop labels
    y = df["IS_ARTIFACT"].copy() # labels

    df_transformed = full_pipeline.fit_transform(X)

    # The arrays passed to models.
    X = df_transformed
    #y = df.IS_ARTIFACT
    assert(X.shape[0] == y.shape[0])
    
    return X, y

# Evaluate model for feature removal
def evaluate_feature_removal_model(model, X, y):
    scores = cross_val_score(model, X, y, scoring="accuracy", cv=3)
    return scores.mean()

# Investigate predictive power of features
# Randomly selects features  

features_and_scores = []

for i in range(1000):
    try:
        forest_cls = RandomForestClassifier(
        n_estimators=100,
        max_depth=7, #4
        random_state=42,
        n_jobs=16,
        class_weight='balanced',
        bootstrap=False #bootstrap true previously
    )

        features_to_keep = select_random_features(3) #5 for ML_df.pkl
        X, y = preprocess_with_features(df, features_to_keep)
        score = evaluate_feature_removal_model(forest_cls, X, y)

        key = tuple(sorted(features_to_keep))
        features_and_scores.append((key, score))
        print(i, key, score)
    except Exception as e:
        print(e)
   

with open('pp9_joined_LR_feats_bsF.pkl', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    # Keeps object type during pickling
    pickle.dump(features_and_scores, f, pickle.HIGHEST_PROTOCOL)   