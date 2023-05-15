This repository stores the codes for the classificiations we have tested for heart murmur detection.


Data
Heart sound recordings of patients are located in data/training_data folder.


Feature Extractions<br />
Window_entropy.m extracts entropy of heart sound signals in the wavelet domain in a rolling window base(window length of 1024).<br />
Window_Monofractal.m compute Hurst exponent in heart sounds from the wavelet domain in a rolling window base(window length of 1024).<br />
Multifractal.m extracts multifractal properties of a heart sounding recording of a patient.


Features
Features extracted above are save in csv files, and csv files are found in case_features for case(murmur present) and control_features for control(murmur absent).


Demonstration
demo_classification.py demonstrates the best model we obatined in our experimentation.


Functions
data_mining_process.py contains functions which are used frequently in our codes in the data mining pipe line.


Experimentation
Train and test classification models with different hyper parameters and different data processing(i.e., balancing controls and cases). Files:
Window_GridSearch_LogisticRegression.ipynb
Window_GridSearch_KNN.ipynb
Window_GridSearch_SVM.ipynb
Window_features_NN.ipynb
