function u = TestOutcome(ytest, ypred)
    tp = 0; tn = 0; fp = 0; fn = 0;
    for i = 1: length(ypred)
        if ytest(i) == 1 && ypred(i) == 1
            tp = tp + 1;
        elseif ytest(i) == 0 && ypred(i) == 1
            fp = fp + 1;
        elseif ytest(i) == 1 && ypred(i) == 0
            fn = fn + 1;
        elseif ytest(i) == 0 && ypred(i) == 0
            tn = tn + 1;
        end
    end
    u = [tp, fp, fn, tn];
end 
