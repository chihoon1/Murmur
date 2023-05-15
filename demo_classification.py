'''
This files demonstrates the best classifier model(SVM) obtained with our approach in wavelet domain
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import MinMaxScaler,StandardScaler
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
from sklearn.metrics import recall_score, precision_score, f1_score
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.utils import class_weight

import warnings
warnings.filterwarnings('ignore')

# make directory for visualization
import os
if not os.path.exists('visualization'):
    os.makedirs('visualization')

# load cases from case features folder
Ent_Case = pd.read_csv('case_features/Window_Entropy_Case.csv', header=None)
Slope_Case = pd.read_csv('case_features/Window_Slope_Case.csv', header=None)
Mf_Case = pd.read_csv('case_features/Window_Mf_Case.csv', header=None)

# load controls from control features folder
Ent_Control = pd.read_csv('control_features/Window_Entropy_Control.csv', header=None)
Slope_Control = pd.read_csv('control_features/Window_Slope_Control.csv', header=None)
Mf_Control = pd.read_csv('control_features/Window_Mf_Control.csv', header=None)

col_names = ["Entropy ", "Slope ", "ID", "Hurst Exponent ",
    "Left Slope ", "Right Slope ","Left Tangent ",
    "Right Tangent", "Broadness", "Left Tangent Point", "Right Tangent Point"]
# combine Ent_Case and Slope_Case and Mfcc_Case
Case = pd.concat([ Ent_Case.T, Slope_Case.T, Mf_Case ], axis=1)
Case.columns = col_names
Control = pd.concat([ Ent_Control.T, Slope_Control.T, Mf_Control ], axis=1)
Control.columns = col_names
Case['Case'] = 1
Control['Case'] = 0

# combine Case and Control data
data = pd.concat([Case, Control], axis=0)

# check if any missing values
data.isna().sum()

# suffle data
data = data.sample(frac=1, random_state= 42).reset_index(drop=True);
# replace NaN values with mean of the column
data = data.fillna(data.mean())

# Using all features in the data
target = data.columns[11]  # target variable name
features = [col_name for col_name in data.columns[1:11]]  # predictor features name
features.pop(features.index('ID'))  # remove patient ID from feature names list
features.pop(features.index('Broadness'))  # not discriminatory according to wilcoxon ranksum test
features.pop(features.index('Right Tangent'))  # not discriminatory according to wilcoxon ranksum test
features

# split data into training and testing data
X = data.loc[:,features]
y = data.loc[:, target]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# Compute the class weights
class_weights = class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(y_train), y=y_train)
class_weights_dict = dict(enumerate(class_weights))

# initialize a model with the optimal hyperparameters found in our experiment in Window_GridSearch_SVM file
svm = SVC(probability=True, C=100, class_weight=class_weights_dict, gamma=0.1, kernel='rbf')
svm.fit(X_train, y_train)

# performance report
# train accuracy
train_pred = svm.predict(X_train)
print(f"Train Accuracy: {accuracy_score(y_train, train_pred)}")
# performance on test data
y_pred = svm.predict(X_test)
cm = confusion_matrix(y_test, y_pred)
cm = pd.DataFrame({"Predicted Negative(Absent)": cm[:, 0], "Predicted Positive(Present)": cm[:, 1]})
cm.index = ["Actual Negative(Absent)", "Actual positive(Present)"]
print(classification_report(y_test, y_pred))
y_prob = svm.predict_proba(X_test)
print(f"AUC: {roc_auc_score(y_test, y_prob[:, 1])}")
print(f"Test Accuracy: {accuracy_score(y_test, y_pred)}")
print(f"Present F1 score: {f1_score(y_test, y_pred)}")
print(f"Absent F1 score: {f1_score(y_test, y_pred,pos_label=0)}")
print(f"Sensitivity: {recall_score(y_test, y_pred)}")
print(f"Speicificity: {recall_score(y_test, y_pred,pos_label=0)}")
print(f"Confusion Matrix:\n{cm}")
