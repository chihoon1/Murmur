function Entropy = info_entropy (input, b)
% In information theory, entropy is a measure of the uncertainty associated 
% with a random variable. In this context, the term usually refers to the 
% Shannon entropy, which quantifies the expected value of the information 
% contained in a message, usually in units such as bits. In this context, a 
% 'message' means a specific realization of the random variable.
% Shannon denoted the entropy H of a discrete random variable X with possible 
% values {x1, ..., xn} as,
%
%     H(X) = E(I(X)).
% 
% E is the expected value,
% I is the information content of X
% 
% I(X) is itself a random variable. If p denotes the probability mass function 
% of X then the entropy can explicitly be written as,
%             n                n                         n
%     H(X) = SUM p(xi)I(xi) = SUM p(xi)logb(1/p(xi)) = -SUM p(xi)logb(p(xi))
%            i=1              i=1                       i=1
%
% Example usage
% -------------
% in = [.25 .25 .25 .25];
% b = 'bit';
% Entropy = info_entropy (in, b)
Entropy = 1;
if b == 'bit'
    info_sz = size (input);
    n = info_sz (1, 2);
    info_sum = sum (input);
    %if info_sum > 1
    %    error ('Sum of probabilities cannot be greater than 1');
    %end
    Entropy = 0;
    for i=1:n
        if isequal(input(1,i),0)
            tmp = 0;
        else
            tmp = (-input(1,i)*log2(input(1,i)));
        end
        Entropy = Entropy + tmp;
    end
        
elseif b == 'nat'
    info_sz = size (input);
    n = info_sz (1, 2);
    info_sum = sum (input);
    if info_sum > 1
        error ('Sum of probabilities cannot be greater than 1');
    end
    Entropy = 0;
    for i=1:n
        if isequal(input(1,i),0)
            tmp = 0;
        else
            tmp = (-input(1,i)*log(input(1,i)));
        end
        Entropy = Entropy + tmp;
    end
        
elseif b == 'dit'
    info_sz = size (input);
    n = info_sz (1, 2);
    info_sum = sum (input);
    if info_sum > 1
        error ('Sum of probabilities cannot be greater than 1');
    end
    Entropy = 0;
    for i=1:n
        if isequal(input(1,i),0)
            tmp = 0;
        else
            tmp = (-input(1,i)*log10(input(1,i)));
        end
        Entropy = Entropy + tmp;
    end
        
end
        