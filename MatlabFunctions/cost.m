function [J grad h th] = cost(theta, xtrain, ytrain, alpha, iter);
    th = theta;
     m = size(xtrain,1);
    for j = 1:iter
        h = sigmoid(xtrain*th);
        J = -(1/m)*sum( ytrain.*log(h) + (1-ytrain).* log(1-h) );
       th = th + ( alpha/length(xtrain) ) * xtrain' * (ytrain-h);
    end
    grad = zeros(size(theta,1),1);
    
    for i = 1:size(grad)
        grad(i) = (1/m)*sum((h-ytrain)'*xtrain(:,i));
    end
end

