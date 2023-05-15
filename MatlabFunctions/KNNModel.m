function [acc_ktr, acc_ts, sensi, speci] = KNNModel(Data, k)  

X_tr = Data.X_tr; Y_tr = Data.Y_tr; 
X_ts = Data.X_ts; Y_ts = Data.Y_ts;

%X_tr1 = zscore(X_tr);% Normalized Predictors
X_tr1 = X_tr;
%X_ts1 = zscore(X_ts);% Normalized Predictors
X_ts1 = X_ts;

    %KnnModel = fitcknn(X_tr1,Y_tr,'NumNeighbors',k, 'Standardize',1);
    KnnModel = fitcknn(X_tr1, Y_tr,'NumNeighbors',k,...
    'Distance','mahalanobis',...
    'Standardize',1);
    
    % KNN validation
    [p3, q3] = predict(KnnModel, X_tr1);
    acc_ktr = (1 - sum(abs(p3 - Y_tr))/length(Y_tr))*100;% correct classification percentage
    
    [p4, q4] = predict(KnnModel, X_ts1);
    acc_ts = (1 - sum(abs(p4 - Y_ts))/length(Y_ts))*100;% correct classification percentage
    
    % sensitivity and specificivity 
    E_k = TestOutcome(Y_ts, p4);
    [sensi,speci] = sesp(E_k(1),E_k(2), E_k(3), E_k(4));
    %[se_k sp_k pre_k ppv_k npv_k ac_k] = sesp(E_k(1),E_k(2), E_k(3), E_k(4));
end 