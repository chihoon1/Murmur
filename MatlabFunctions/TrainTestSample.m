
function [X_train, y_train, X_test, y_test] = TrainTestSample(X_control, X_case, tr_per)
X_c2 = X_case; 
X_n2 = X_control;


y_c = repelem(1, size(X_c2,1)); y_n = repelem(0,size(X_n2,1)); % asign labels case and control as 1 and 0

A = min(length(y_c), length(y_n));% select number of samples for fitting a classification model

ids = randperm(size(y_n,2), A); 
X_n3 = X_n2(ids,:); y_n3 = y_n(ids);

X = cat(1, X_n3, X_c2); y = cat(2, y_c, y_n3);

cv = cvpartition(size(X,1),'HoldOut',tr_per);
idx = cv.test;

X_train = X(idx,:); X_test = X(~idx, :);
y_train = y(idx)'; y_test = y(~idx)';


