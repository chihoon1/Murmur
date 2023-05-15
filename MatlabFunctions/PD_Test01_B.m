close all; clear all; clc
addpath('/Users/dixon/Documents/TAMU/DemosNew/')
addpath('/Users/dixon/Documents/TAMU/Sample Data/White blood Cancer data/MatlabFunctions/')
%% INITIALIZING PARAMETERS 
J = 13;              % power of two used for selecting data points
n = 2^J;             % number of data points selected from data files 
L = 1;               % wavelet decomposition level
k1 = L+1; k2 = J - 1;  % range of scale used to compute wavelet spectra 
ismean = 0;          % measure used to compute wavelet spectra; options 0- mean/ 1 - median
isplot = 1;          % plot wavelet spectra 0 - No/ 1 - yes
Lt = 18; Rt = 19;    % column index coressponding to left and right foot in data file
filt = [sqrt(2)/2 sqrt(2)/2];   % wavelet filter 
     
%% 
% ##############################  Control Data Processing ####################
% (1) Read Row Data files into Matlab
dirName = sprintf('/Users/dixon/Documents/TAMU/Parkinson Disease/Data/Control/');             %# folder path
files = dir( fullfile(dirName,'*.txt') );   %# list all *.xyz files
files = {files.name}';                      %'# file names
nfi = numel(files);

% (2) Compute wavelet spectra and then derive slope
L_Slope = []; R_Slope = [];
for i = 1 : nfi
    fname = fullfile(dirName,files{i});     %# full path to file
    data = readmatrix(fname);
    
    if size(data,1) >= n
        Ls = []; Rs = [];
        for j = 2 : 9
            % lf - left foot 
            lf = zscore(data(1:n+1, j ));
            lf = cumsum(lf);
            [slope_lf] = waveletspectra(lf, L, filt, k1, k2, ismean, isplot);
            Ls = [Ls slope_lf];

            % rf - left foot 
            rf = zscore(data(1:n+1,j + 8));
            rf = cumsum(rf);
            [slope_rf] = waveletspectra(rf, L, filt, k1, k2, ismean, isplot);
            Rs = [Rs slope_rf];
        end 
        
    L_Slope = [L_Slope; median(Ls) ]; R_Slope = [R_Slope; median(Rs) ];
    end 
end 

Control_Slope = [L_Slope R_Slope];

%(3) Compute level-wise covariance ( Wavelet Flux )
Control_CV = [];
for i = 1 : nfi
    fname = fullfile(dirName,files{i});     %# full path to file
    data = readmatrix(fname);
       
    if size(data,1) >= n
        CV = [];
        for p = 2:9
            %lf - left foot 
            lf = data(1:n, p); wc_lf = dwtr(lf, filt, J - L);
            %rf - left foot
            rf = data(1:n, p+8); wc_rf = dwtr(rf, filt, J - L);

            CV_LR = [];
            for j = L : J - 1
                help_lf = wc_lf( round(2^(j)+1) : round(2^(j+1)) );  
                help_rf = wc_rf( round(2^(j)+1) : round(2^(j+1)) ); 

                st_lf = std(help_lf); st_rf = std(help_rf);
                
                cv = help_lf*help_rf'/ (st_lf * st_rf);
                
                CV_LR = [CV_LR cv];
            end
            CV = [ CV; CV_LR];
        end 
    Control_CV = [Control_CV; median(CV,1) ];
    
    end 
end 

% %%%% ###########################  Case Data Processing #################
% %
%% 
% ##############################  Control Data Processing ####################
% (1) Read Row Data files into Matlab
dirName = sprintf('/Users/dixon/Documents/TAMU/Parkinson Disease/Data/Case/');             %# folder path
files = dir( fullfile(dirName,'*.txt') );   %# list all *.xyz files
files = {files.name}';                      %'# file names
nfi = numel(files);

% (2) Compute wavelet spectra and then derive slope
L_Slope = []; R_Slope = [];
for i = 1 : nfi
    fname = fullfile(dirName,files{i});     %# full path to file
    data = readmatrix(fname);
    
    if size(data,1) >= n
        Ls = []; Rs = [];
        for j = 2 : 9
            % lf - left foot 
            lf = zscore(data(1:n+1, j ));
            lf = cumsum(lf);
            [slope_lf] = waveletspectra(lf, L, filt, k1, k2, ismean, isplot);
            Ls = [Ls slope_lf];

            % rf - left foot 
            rf = zscore(data(1:n+1,j + 8));
            rf = cumsum(rf);
            [slope_rf] = waveletspectra(rf, L, filt, k1, k2, ismean, isplot);
            Rs = [Rs slope_rf];
        end 
        
    L_Slope = [L_Slope; median(Ls) ]; R_Slope = [R_Slope; median(Rs) ];
    end 
end 

Case_Slope = [L_Slope R_Slope];

%(3) Compute level-wise covariance ( Wavelet Flux )
Case_CV = [];
for i = 1 : nfi
    fname = fullfile(dirName,files{i});     %# full path to file
    data = readmatrix(fname);
       
    if size(data,1) >= n
        CV = [];
        for p = 2:9
            %lf - left foot 
            lf = data(1:n, p); wc_lf = dwtr(lf, filt, J - L);
            %rf - left foot
            rf = data(1:n, p+8); wc_rf = dwtr(rf, filt, J - L);

            CV_LR = [];
            for j = L : J - 1
                help_lf = wc_lf( round(2^(j)+1) : round(2^(j+1)) );  
                help_rf = wc_rf( round(2^(j)+1) : round(2^(j+1)) ); 

                st_lf = std(help_lf); st_rf = std(help_rf);
                
                cv = help_lf*help_rf'/ (st_lf * st_rf);
                
                CV_LR = [CV_LR cv];
            end
            CV = [ CV; CV_LR];
        end 
    Case_CV = [Case_CV; median(CV,1) ];
    
    end 
end 
%
%% Contruct feature matrices and data preparation for classification

X_case = cat(2, Case_Slope, Case_CV); X_control = cat(2, Control_Slope, Control_CV);
Y_case = repelem(1, size(X_case,1));       Y_control = repelem(1, size(X_control,1));

F = (median(X_case,1) - median(X_control,1)).^2 ./ (var(X_case,1) + var(X_control,1));

[q, r ]= maxk(F, 5);

X_case = X_case(:,r); X_control = X_control(:,r);


%% Classification 

% number of repititions 
nsample = 1000;
% divide dataset into training and testing data 
tr_per = .80; % percentrage of samples used to train classifiers 

A = zeros(nsample,3,4);

for i = 1 : nsample
[X_tr, Y_tr, X_ts, Y_ts] = TrainTestSample(X_control, X_case, tr_per);

% training data set
Data.X_tr = X_tr; Data.Y_tr = Y_tr;

% testing dataset 
Data.X_ts = X_ts; Data.Y_ts = Y_ts; 

% (1) Logistic regression 
p = 0.5; % threshold
[acc_tr, acc_ts, sensi, speci] = LogisticModel(Data, p);

A(i,1,:) = [acc_tr, acc_ts, sensi, speci];

% (2) Support Vector Machine

[acc_str, acc_ts, sensi, speci] = SVMMOdel(Data);
A(i,2,:) = [acc_str, acc_ts, sensi, speci];

% (3) k-Nearest neightbor 
k = 5;
[acc_ktr, acc_ts, sensi, speci] = KNNModel(Data, k);
A(i,3,:) = [acc_ktr, acc_ts, sensi, speci];

end 
k = zeros(3,4);
k(:,:) = mean(A,1) 


% % 
% %%
%X_case = X_case(:,[1 2 3 4 5]); X_control = X_control(:,[1 2 3 4 5]);
y_ca = repelem(1, size(X_case,1)); y_co = repelem(0,size(X_control,1));
 
X = zscore([X_case; X_control]);
Y = [y_ca'; y_co'];
 
a = randperm(length(Y));
 
X = X(a,:); Y = Y(a); 
T = table( X(:,1), X(:,2),X(:,3), X(:,4), X(:,5), Y);
 
rng("default") % For reproducibility of the partition
c1 = cvpartition(T.Y,"Holdout",0.33);
trainingIndices1 = training(c1);
validationIndices1 = test(c1);

tblTrain1 = T(trainingIndices1,:);
tblValidation1 = T(validationIndices1,:);

Mdl = fitcnet(tblTrain1,"Y", "ValidationData",tblValidation1, ...
    "Verbose",1, "LayerSizes",[5]);
%"LayerSizes",[ 5 2], , "LayerSizes",[ 10 7]


figure(2)
subplot(311)
[labels,Scores] = predict(Mdl,tblTrain1);
confusionchart(tblTrain1.Y,labels)
subplot(312)
[labels,Scores] = predict(Mdl,tblValidation1);
confusionchart(tblValidation1.Y,labels)


error = loss(Mdl,tblTrain1,"Y");
Train_accuracy = (1-error)*100

error = loss(Mdl,tblValidation1,"Y");
Test_accuracy = (1-error)*100


subplot(313)
iteration = Mdl.TrainingHistory.Iteration;
trainLosses = Mdl.TrainingHistory.TrainingLoss;
valLosses = Mdl.TrainingHistory.ValidationLoss;

plot(iteration,trainLosses)
hold on 
plot(iteration,valLosses)
legend(["Training","Validation"])
xlabel("Iteration")
ylabel("Cross-Entropy Loss")

hold off
% %% Decision Tree
% 
% rng('default') 
% tallrng('default')
% Mdl = fitctree(X(1:200,:),Y(1:200));
% 
% cvMdl = crossval(Mdl);
% L = kfoldLoss(cvMdl)
% 
% figure(3)
% subplot(121)
% [labls,Scores] = predict(Mdl,X(1:200,:));
% confusionchart(Y(1:200),labls)
% 
% 
% subplot(122)
% [labls,Scores] = predict(Mdl,X(200:end,:));
% confusionchart(Y(200:end),labls)