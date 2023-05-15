function [se sp pre ppv npv ac] = sesp(tp, fp, fn, tn)
%D-dimer as a test for acute PE (Goldhaber et al, 1993)
% [s1, s2, p1, p2, p3, a, yi] = sesp(42,96,3,32)
%
n = tp+tn+fn+fp; %total sample size
np = tp + fp; %total positive
nn = tn + fn; %total negative
nd = tp + fn; %total with disease
nc = tn + fp; %total control (without disease)
%--------------
se = tp/nd; %tp/(tp + fn):::sensitivity
sp = tn/nc; %tn/(tn + fp):::specificity
pre = nd/n; %(tp + fn)/(tp+tn+fn+fp):::prevalence
%only in for case when sample is random from the
% population of interest. Otherwise, the prevalence
% needed for calculating PPV and NPV is an input value
ppv = tp/np; %tp/(tp + fp):::positive predictive value
npv = tn/nn; %tn/(tn+fn):::negative predictive value
lrp = se/(1-sp); %:::likelihood ratio positive
lrn = (1-se)/sp; %:::likelihood ratio negative
ac = (tp+tn)/n; %:::accuracy
yi = (se + sp - 1)/sqrt(2); %:::youden indexsesp(42,96,3,32)
%---------------
%disp(' Se Sp Pre PPV NPV LRP LRN Ac Yi')
%disp([se, sp, pre, ppv, npv, lrp, lrn, ac, yi])