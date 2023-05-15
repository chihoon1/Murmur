function [acc_str, acc_ts, sensi, speci] = SVMMOdel(Data)
X_tr = Data.X_tr; Y_tr = Data.Y_tr; 
X_ts = Data.X_ts; Y_ts = Data.Y_ts;

    %X_tr1 = zscore(X_tr);% Normalized Predictors
    X_tr1 = X_tr;
    %SVMModel = fitcsvm(X_tr1,Y_tr,'KernelFunction','linear','Standardize',false);
    
    SVMModel = fitcsvm(X_tr1,Y_tr,'KernelFunction','linear','Standardize',false, 'Solver', 'SMO',...
        'KernelScale',0.0014051  , 'BoxConstraint', 0.0010045 );
  
    % SVM validation   
    [p1, q1] = predict(SVMModel, X_tr1);
    acc_str = (1 - sum(abs(p1 - Y_tr))/length(Y_tr))*100;% correct classification percentage
     
    %X_ts1 = zscore(X_ts);% Normalized Predictors
    X_ts1 = X_ts;
    [p2, q2] = predict(SVMModel, X_ts1);
    
    acc_ts = (1 - sum(abs(p2 - Y_ts))/length(Y_ts))*100;% correct classification percentage
    
     % sensitivity and specificivity 
     E_s = TestOutcome(Y_ts, p2);
     %[se_s,sp_s] = sesp(E_s(1),E_s(2), E_s(3), E_s(4));
     [sensi, speci, pre_s, ppv_s, npv_s, ac_s] = sesp(E_s(1),E_s(2), E_s(3), E_s(4));
end 