#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import os.path
#import sklearn as sk
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

####################

df = pd.read_pickle('ML_df.pkl')

####################

X = df.drop(["IS_ARTIFACT"], axis=1) # drop labels
y = df["IS_ARTIFACT"].copy() # labels

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create testing and training set
#train_set, test_set = train_test_split(df, test_size=0.2, random_state=42)

# Seperate predictors and labels in training set to transform the training attributes without affecting the predictors
#artifact_train = train_set.drop(["IS_ARTIFACT"], axis=1) # drop labels for training set 
#artifact_train_labels = train_set["IS_ARTIFACT"].copy() # IS_ARTIFACT is the label (boolean)

########################
####  Data Cleaning ####
########################

# 1. Handle numerical attributes
imputer = SimpleImputer(strategy='median')
cat_attribs = ['REF','ALT'] 
# Need to remove text attribute bc median can only be calculated on numerical attributes
X_train_num = X_train.drop(cat_attribs, axis=1) 
# Now we can fit the imputer instance to the training data using the fit() method
imputer.fit(X_train_num)
# Now we can use the "trained" imputer to tranform the training set by replacing missing values by the learned median
# Transform the training set:
X = imputer.transform(X_train_num)

# 2. Handle text/categorical attributes
X_train_cat = X_train[cat_attribs]
# One-hot code categorical attributes
cat_encoder = OneHotEncoder()
X_train_cat_1hot = cat_encoder.fit_transform(X_train_cat)

###################################
#### Put together the pipeline ####
###################################

# Pipeline for numerical attributes
num_pipeline =  Pipeline([
    ('imputer', SimpleImputer(strategy="median")), # fill in missing values with median values
    ('std_scaler', StandardScaler()),              # feature scaling
])

artifact_num_tr = num_pipeline.fit_transform(X_train_num)

# Complete the pipeline with combined categorical attributes
num_attribs = list(X_train_num)
# cat_attribs = ['REF','ALT'] defined previously in step 1 already

full_pipeline = ColumnTransformer([
            ("num", num_pipeline, num_attribs),
            ("cat", OneHotEncoder(), cat_attribs),
])

X_train_prepared = full_pipeline.fit_transform(X_train)
# Transform test set, NOT fit_transform bc we dont want to fit the test set
X_test_prepared = full_pipeline.transform(X_test)


# lazypredict.log has a botleneck model that takes long to build

# lazypredict1.log
# same list as lazypredict default gives but 
# removed ('NuSVC', sklearn.svm._classes.NuSVC)

# lazypredict2.log
# remove ('SVC', sklearn.svm._classes.SVC) bc another bottleneck

clf_list = [('AdaBoostClassifier', sklearn.ensemble._weight_boosting.AdaBoostClassifier),
 ('BaggingClassifier', sklearn.ensemble._bagging.BaggingClassifier),
 ('BernoulliNB', sklearn.naive_bayes.BernoulliNB),
 ('CalibratedClassifierCV', sklearn.calibration.CalibratedClassifierCV),
 ('CategoricalNB', sklearn.naive_bayes.CategoricalNB),
 ('DecisionTreeClassifier', sklearn.tree._classes.DecisionTreeClassifier),
 ('DummyClassifier', sklearn.dummy.DummyClassifier),
 ('ExtraTreeClassifier', sklearn.tree._classes.ExtraTreeClassifier),
 ('ExtraTreesClassifier', sklearn.ensemble._forest.ExtraTreesClassifier),
 ('GaussianNB', sklearn.naive_bayes.GaussianNB),
 ('KNeighborsClassifier',
  sklearn.neighbors._classification.KNeighborsClassifier),
 ('LabelPropagation',
  sklearn.semi_supervised._label_propagation.LabelPropagation),
 ('LabelSpreading', sklearn.semi_supervised._label_propagation.LabelSpreading),
 ('LinearDiscriminantAnalysis',
  sklearn.discriminant_analysis.LinearDiscriminantAnalysis),
 ('LinearSVC', sklearn.svm._classes.LinearSVC),
 ('LogisticRegression', sklearn.linear_model._logistic.LogisticRegression),
 ('NearestCentroid', sklearn.neighbors._nearest_centroid.NearestCentroid),
 ('PassiveAggressiveClassifier',
  sklearn.linear_model._passive_aggressive.PassiveAggressiveClassifier),
 ('Perceptron', sklearn.linear_model._perceptron.Perceptron),
 ('QuadraticDiscriminantAnalysis',
  sklearn.discriminant_analysis.QuadraticDiscriminantAnalysis),
 ('RandomForestClassifier', sklearn.ensemble._forest.RandomForestClassifier),
 ('RidgeClassifier', sklearn.linear_model._ridge.RidgeClassifier),
 ('RidgeClassifierCV', sklearn.linear_model._ridge.RidgeClassifierCV),
 ('SGDClassifier', sklearn.linear_model._stochastic_gradient.SGDClassifier),
 ('StackingClassifier', sklearn.ensemble._stacking.StackingClassifier),
 ('XGBClassifier', xgboost.XGBClassifier),
 ('LGBMClassifier', lightgbm.LGBMClassifier)]

#lazypredict.Supervised.CLASSIFIERS same as clf_list except XGB and LGBM
#('XGBClassifier', xgboost.sklearn.XGBClassifier)
#('LGBMClassifier', lightgbm.sklearn.LGBMClassifier)

clf = LazyClassifier(verbose=2,ignore_warnings=True, custom_metric=None, classifiers=clf_list)
#clf = LazyClassifier(verbose=0,ignore_warnings=True, custom_metric=None)

models,predictions = clf.fit(X_train_prepared, X_test_prepared, y_train, y_test)


models.to_csv('models2.txt', sep='\t')
predictions.to_csv('predictions2.txt', sep='\t')


#models.to_csv('models.txt', sep='\t')
#predictions.to_csv('predictions.txt', sep='\t')