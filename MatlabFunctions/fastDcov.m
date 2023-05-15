function covsq = fastDcov(x,y)

n = length(x);
% n = 4;
% x = [5; 7; 9; 8];
% y = [1; 2; 3; 4];
[x, Index] = sort(x);
y = y(Index);

si = cumsum(x);
s = si(n);
a_x = (-(n-2):2:n).' .*x + (s - 2*si);


v = [x y x.*y];
nw = size(v, 2);
idx = zeros(n, 2);
idx(:,1) = 1:n;

%idx_iv = 2:(nw + 1);
%iv = zeros(n, nw + 1);
iv1 =  zeros(n, 1);
iv2 =  zeros(n, 1);
iv3 =  zeros(n, 1);
iv4 =  zeros(n, 1);

i = 1; r = 1; s = 2;
while i < n;
    gap = 2*i;
    k = 0;
    idx_r = idx(:,r);
    csumv = [zeros(1, nw); cumsum(v(idx_r,:))];
    
    for j = 1:gap:n;
        st1 = j; e1 = min(st1 + i - 1,n);
        st2 = j + i; e2 = min(st2 + i - 1,n);
        
        while (st1 <= e1) && (st2 <= e2);
            k = k +1;
            idx1 = idx_r(st1);
            idx2 = idx_r(st2);
            
            if y(idx1) >= y(idx2);
                idx(k,s) = idx1;
                st1 = st1 + 1;
            else
                idx(k,s) = idx2;
                st2 = st2 + 1;
                %iv(idx2, 1) = iv(idx2, 1) + e1 -st1 +1;
                %iv(idx2, idx_iv) = iv(idx2, idx_iv) + (csumv(e1+1, :) - csumv(st1, :));
                %iv(idx2, 2) = iv(idx2, 2) + (csumv(e1+1, 1) - csumv(st1, 1));
                %iv(idx2, 3) = iv(idx2, 3) + (csumv(e1+1, 2) - csumv(st1, 2));
                %iv(idx2, 4) = iv(idx2, 4) + (csumv(e1+1, 3) - csumv(st1, 3));
                
                iv1(idx2, 1) = iv1(idx2) + e1 -st1 +1;
                %iv(idx2, idx_iv) = iv(idx2, idx_iv) + (csumv(e1+1, :) - csumv(st1, :));
                iv2(idx2) = iv2(idx2) + (csumv(e1+1, 1) - csumv(st1, 1));
                iv3(idx2) = iv3(idx2) + (csumv(e1+1, 2) - csumv(st1, 2));
                iv4(idx2) = iv4(idx2) + (csumv(e1+1, 3) - csumv(st1, 3));
            end;
          
        end;   
        
        if st1 <= e1;
            kf = k + e1 - st1 + 1;
            idx((k+1):kf, s) = idx_r(st1:e1,:);
            k = kf;
        elseif st2 <=e2;
            kf = k + e2 - st2 + 1;
            idx((k+1):kf, s) = idx_r(st2:e2,:);
            k = kf;
        end;
            
        
    end;
    
    i = gap;
    r = 3- r; s = 3 -s;
end;

% Caculate b_y
ySorted = y(idx(n:-1:1, r));
si = cumsum(ySorted);
s = si(n);
b_y = zeros(n, 1);
b_y(idx(n:-1:1, r)) = (-(n-2):2:n).' .*ySorted + (s - 2*si);

% caculate D
covterm =  n*(x - mean(x)).' * (y - mean(y));
% c1 = iv(:, 1).' * v(:, 3);
% c2 = sum(iv(:,4));
% c3 = iv(:, 2).' * y;
% c4 = iv(:, 3).' * x;
c1 = iv1.' * v(:, 3);
c2 = sum(iv4);
c3 = iv2.' * y;
c4 = iv3.' * x;

d = 4*((c1 + c2) - (c3 + c4)) - 2*covterm;

nsq = n*n;
ncb = nsq*n;
nq = ncb*n;
term1 = d / nsq;
term2 = 2* (a_x.' * b_y) / ncb;
term3 = sum(a_x) * sum(b_y) / nq;

covsq = (term1 + term3) - term2;

