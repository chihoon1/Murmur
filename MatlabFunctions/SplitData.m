function Y = SplitData(y, L)
N = length(y);  %this varies for different data X
%L = 2048; % this is fixed
X = y;%rand(N,1);
nSplits = ceil(N/L); % we need nSplit segments
sz=L*ones(nSplits,1); % assume all are full, sizes for the segments
sz(end)=N-(nSplits-1)*L; %Fix size of last one, in case less than L are left over
Y=mat2cell(X,sz); %Split it