import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score, recall_score, precision_score, f1_score, confusion_matrix, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from sklearn.model_selection import validation_curve, train_test_split, GridSearchCV
from sklearn.preprocessing import MinMaxScaler



warnings.filterwarnings("ignore")


def MinMaxSclProcess(data, train_data, test_data):
    # Parameters
    # data: dataset of interest with only selected features from the original data set without the target variable
    # train_data: splitted data(param) for training. test_data splitted data(param) for testing
    # Output: MinMax scaled train data and test data
    scaler = MinMaxScaler()
    scaler.fit(data)
    scaled_X_train = pd.DataFrame(scaler.transform(train_data))
    scaled_X_train.columns = data.columns
    scaled_X_test = pd.DataFrame(scaler.transform(test_data))
    scaled_X_test.columns = data.columns
    return scaled_X_train, scaled_X_test


def extract_data_by_classes(data, target_name, predictor_features):
    # param: data(pandas DataFrame) is the dataset
    # param: target_name(string) is the name of the target variable
    # param: predictor_features(list) is the list of columns(variables) that will predict the target variable
    # return: list of dataframes each of which is filtered by a class in the set of classes in the target variable
    class_labels = data[target_name].unique()
    dfs_by_class = []
    columns = predictor_features + [target_name]
    for cls in class_labels:
        dfs_by_class.append(data[data[target_name] == cls].loc[:, columns])
    return dfs_by_class


def sampling_data(dfs_by_class):
    # param: dfs_by_class is a list of pandas dataframes
    # each of which is filtered by a class in the set of classes in the target variable
    # return: one merged data consisting equal number of randomly sampled data points from each class
    min_num_rows = float('inf')
    for df in dfs_by_class:
        num_rows = df.shape[0]
        # finds the minimum number of data points in each dataframe
        if num_rows < min_num_rows:
            min_num_rows = num_rows
    # random sample the data such that each dataframe filtered by a class contains the equal number of data points
    # this prevents training bias toward the overrepresented class
    random_seed = np.random.randint(1,1000)
    sampled_merged_data = dfs_by_class[0].sample(n=min_num_rows, random_state=random_seed)
    for i in range(1,len(dfs_by_class)):
        sampled_class_i_df = dfs_by_class[i].sample(n=min_num_rows, random_state=random_seed)
        sampled_merged_data = sampled_merged_data.append(sampled_class_i_df)
    return sampled_merged_data.sort_index()


def drop_missing_values(data, column_name, missing_value=np.nan):
    # param: data(pandas DataFrame) is the dataset
    # param: column_name(str) is a name of column where missing values are
    # param: missing_value(int, float, str, etc) is the value of missing value. Default numpy.nan
    # return new dataframe where rows with missing value in the column_name are dropped
    return data[data[column_name]!=missing_value]


def repeat_sampling_and_training(model_function, model_f_params, data,
                                  target_name, predictor_features, num_repeat=1000, doMinMaxScaling=False):
    # param: model_function is a function object that creates a model
    # param: model_f_params is a list specifying the parameters of model_function (the order of elements in this list matters)
    # param: data(pandas DataFrame) is the dataset
    # param: target_name(string) is the name of the target variable
    # param: predictor_features(list) is the list of columns(variables) that will predict the target variable
    # param: num_repeat(int) is number of times to repeat sampling data and training model from it
    # param: doMinMaxScaling(Boolean) indicates whether to apply MinMaxScale to data
    # return: dictionary containing performance metrics for each iteration
    # and confusion matrix of the model with the best test accuracy
    performance_measures_dict = {
        'training_accuracy': np.zeros(num_repeat),
        'testing_accuracy': np.zeros(num_repeat),
        'sensitivity': np.zeros(num_repeat),
        'specificity': np.zeros(num_repeat),
        'case_precision': np.zeros(num_repeat),
        'control_precision': np.zeros(num_repeat),
        'case_F1_score': np.zeros(num_repeat),
        'control_F1_score': np.zeros(num_repeat),
        'AUC_score': np.zeros(num_repeat)
    }
    max_test_acc, max_cm = 0, None
    for i in range(num_repeat):
        # random sampling data such that the dataset contains the equal number of rows from each class
        dfs_lst = extract_data_by_classes(data, target_name, predictor_features)
        sampled_data = sampling_data(dfs_lst)
        # split into training and testing data

        X = sampled_data.drop(labels=target_name, axis=1)
        y = sampled_data[target_name]
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

        # scaling data
        if doMinMaxScaling:
            X_train, X_test = MinMaxSclProcess(pd.concat([X_train, X_test], axis=0), X_train, X_test)

        # fit model to training data
        model = model_function(model_f_params)
        model.fit(X_train, y_train)

        # predict y values for X_test
        train_pred = model.predict(X_train)
        test_pred = model.predict(X_test)

        # store performance metrics for each iteration
        test_acc = accuracy_score(y_test, test_pred)
        performance_measures_dict['training_accuracy'][i] = accuracy_score(y_train, train_pred)
        performance_measures_dict['testing_accuracy'][i] = test_acc
        performance_measures_dict['sensitivity'][i] = recall_score(y_test, test_pred, pos_label=1)
        performance_measures_dict['specificity'][i] = recall_score(y_test, test_pred, pos_label=0)
        performance_measures_dict['case_precision'][i] = precision_score(y_test, test_pred, pos_label=1)
        performance_measures_dict['control_precision'][i] = precision_score(y_test, test_pred, pos_label=0)
        performance_measures_dict['case_F1_score'][i] = f1_score(y_test, test_pred, pos_label=1)
        performance_measures_dict['control_F1_score'][i] = f1_score(y_test, test_pred, pos_label=0)
        y_prob = model.predict_proba(X_test)
        performance_measures_dict['AUC_score'][i] = roc_auc_score(y_test, y_prob[:, 1])
        if test_acc > max_test_acc:
            max_test_acc = test_acc
            max_cm = confusion_matrix(y_test, test_pred)
    max_cm = pd.DataFrame({"Predicted Negative(Absent)": max_cm[:, 0], "Predicted Positive(Present)": max_cm[:, 1]})
    max_cm.index = ["Actual Negative(Absent)", "Actual positive(Present)"]
    return performance_measures_dict, max_cm


def print_performance_metrics(performance_measures_dict):
    # param: performance_measures_dict(dict) contains performance metrics for each iteration
    print(f"Mean of Training Accuracy: {performance_measures_dict['training_accuracy'].mean()}")
    print(f"Mean of Testing Accuracy: {performance_measures_dict['testing_accuracy'].mean()}")
    print(f"Standard deviation of Testing Accuracy: {performance_measures_dict['testing_accuracy'].std()}")
    print(f"Mean of Sensitivity: {performance_measures_dict['sensitivity'].mean()}")
    print(f"Standard deviation of Sensitivity: {performance_measures_dict['sensitivity'].std()}")
    print(f"Mean of Speicificity: {performance_measures_dict['specificity'].mean()}")
    print(f"Standard deviation of Speicificity: {performance_measures_dict['specificity'].std()}")
    print(f"Mean of Case Precision: {performance_measures_dict['case_precision'].mean()}")
    print(f"Standard deviation of Case Precision: {performance_measures_dict['case_precision'].std()}")
    print(f"Mean of Control Precision: {performance_measures_dict['control_precision'].mean()}")
    print(f"Standard deviation of Control Precision: {performance_measures_dict['control_precision'].std()}")
    print(f"Mean of Case F1: {performance_measures_dict['case_F1_score'].mean()}")
    print(f"Standard deviation of Case F1: {performance_measures_dict['case_F1_score'].std()}")
    print(f"Mean of Control F1: {performance_measures_dict['control_F1_score'].mean()}")
    print(f"Standard deviation of Control F1: {performance_measures_dict['control_F1_score'].std()}")
    print(f"Mean of AUC: {performance_measures_dict['AUC_score'].mean()}")
    print(f"Standard deviation of AUC: {performance_measures_dict['AUC_score'].std()}")
