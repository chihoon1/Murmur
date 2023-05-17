# Wavelet-based Approach for Detecting Heart Murmurs


In this project, we propose a novel method based on self-similiarity property of heart sound signls for a diagnosis of heart murmurs.
1. Monofractal self-similiarty property is obtained by computing wavelet spectra in wavelet domain.
2. Moultifracl self-similiarty properties are found based on multifractal spectra in wavelet domain.
3. A set of localized discriminatory features are constructed and selected via a rolling window method to build classifiers; window size is 1024.
4. Heart murmurs detection performance is evaluated by using four classifiers: Logistic Regression, k-nearrest neighbor, Support Vector Machine, and Neural Network.


### Dataset
1. The dataset consists of heart sounds recorded from 1568 patients (305 (19.5%) cases and 1144 (73%) controls) in four ascultation locations and saved as phonocardiogram (PCG) signals. The dataset contrains three groups: murmur present(case), murmur absent(control), and inconclusive.
2. More information about the dataset is given in the paper, Jorge Oliveira et al, The circor digiscope dataset: From murmur detection
to murmur classification. IEEE Journal of Biomedical and Health Informatics, 26(6):2524–2535, 2022.
3. This project considered heart sound signals of patients with heart murmurs and without it.
4. Heart sound recordings of patients are located in data/training_data folder.


### Codes<br />
#### Matlab codes: Feature Extractions
1. Window_entropy.m extracts entropy of heart sound signals in the wavelet domain in a rolling window base(window length of 1024).<br />
2. Window_Monofractal.m compute Hurst exponent in heart sounds from the wavelet domain in a rolling window base(window length of 1024).<br />
3. Multifractal.m extracts multifractal properties of a heart sounding recording of a patient.<br />
4. Extracted features are saved in csv files, and csv files can be found in case_features for case(murmur present) and control_features for control(murmur absent).


#### Python codes:
1. demo_classification.py demonstrates the best model we obatined in our experimentation.<br />
2. data_mining_process.py contains functions which are used frequently in our codes in the data mining pipe line.<br />
3. Train and test classification models with different hyper parameters and different data processing(i.e., balancing controls and cases) in the following files:<br />
        &nbsp Window_GridSearch_LogisticRegression.ipynb<br />
        &nbsp;Window_GridSearch_KNN.ipynb<br />
        &nbsp;Window_GridSearch_SVM.ipynb<br />
        &nbsp;Window_features_NN.ipynb<br />
