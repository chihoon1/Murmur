close all; clear all; clc
addpath('./MatlabFunctions/');

lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

dirName = sprintf('./data/training_data');             %# folder path
files = dir( fullfile(dirName,'*.wav') );   %# list all *.xyz files
files = {files.name}';                      %'# file names
nfi = numel(files);

filt= [-0.075765714789341  -0.029635527645954   0.497618667632458 ...
          0.803738751805216   0.297857795605542  -0.099219543576935 ...
        -0.012603967262261   0.032223100604071];
    
L = 1; k = 10;
l = 2^k; ismean = 2; isplot =0;
 
H = zeros(nfi, 200); Loc = cell(nfi, 1);
Std_H  = zeros(nfi, 2);

LL = zeros(nfi, 2);  W = [];
for i  = 1 :nfi
    % extract patient ID
    newStr = split(files{i},[" ","_","."]);
    Patient_ID = str2num(newStr{1});
    Murmer_Location = newStr{2};
    
    % Read audio file
    fname = fullfile(dirName,files{i});
    [data, fs] = audioread(fname); 
    
    % standardize data
    data = (data - mean(data))/ max(abs(data));
    
    % store murmur location
    LL(i,:) = [Patient_ID length(data)];
    
    % select sound signal of size power of 2
    N=length(data);  LN= log2(N); NN=floor(LN);
    y= data(1:2^NN);
    
    % split data into windows of size 1024
    yi = SplitData(y, l); W = [W, length(yi)];
    
    % store patient ID in the first column of H
    H(i,1) = Patient_ID;
    
    for j = 1: length(yi)   
        % Estimate Hurst exponent
        [h, levels, log2spec ] = waveletspectra(yi{j}, 1, filt, k-4, k-1, ismean, isplot);
        H(i,j+1) = h;
    end
    Loc{i} =  Murmer_Location;
end  

D = [];
for i = 1 :size(H,1)
    a = H(i, 2:end);
    
    D = [D median(a(find(a)))] ;
    
    a(find(~a)) = mean(a(find(a)));
    
    H(i,2:end)= a;
end 

t = readtable('training_data.csv');

a = t(:,1); b = t(:,8); b = string(table2array(b));

control_ID = find(b == 'Absent'); case_ID = find(b == 'Present'); unknown_ID = find(b == 'Unknown');

a1 = table2array(a); 

co_ID = a1(control_ID); ca_ID = a1(case_ID); un_ID = a1(unknown_ID);

CO =[]; CA = []; Un = [];
Class = zeros(size(LL,1),1); L_new = [LL Class]; 
for i =  1:size(LL,1)
    if ismember(LL(i,1), ca_ID)
        CA = [CA i];
        L_new(i,3) = 1;
    elseif  ismember(LL(i,1), co_ID)
        CO = [CO, i];
    else
        Un = [Un i];
    end                
end

% plotting density plots
fig = figure(3);
fig.Position = [15 10 1800 1500];

X_control = D(CO); X_case = D(CA); X_un = D(Un);

for i = 1:1
    subplot(1,1,i)
    [f1,xi1] = ksdensity(X_case); 
    plot(xi1,f1, 'r-','linewidth', 2);

    hold on
    [f2,xi2] = ksdensity(X_control); 
    plot(xi2,f2,'b--','linewidth', 2);
 
    %xlim([-4, 0])
    legend(["Cases", " Controls"], "fontsize", 15)
    %title(sprintf("Feature %d (p-value = %.2f) ",i, z(3,i)),'fontweight','bold','fontsize',15 )
    ylabel("Probability",'fontsize',14 )
    xlabel('Slope','fontsize',14 )
    
    %a = get(gca,'XTickLabel');
    set(gca,"FontSize",15)

    grid on
    hold off
end

% save into csv files
%saveas(fig,'Figures/Window_Slope.png')

Window_H.Case = X_case; Window_H.Control = X_control; Window_H.Un = X_un;
%save('Window_H.mat','Window_H')

writematrix( X_case,'./case_features/Window_Slope_Case.csv');
writematrix( X_control,'./control_features/Window_Slope_Control.csv');
% writematrix( X_un,'Window_Slope_Un.csv');

pt = unique(H(:,1)); N = length(pt);  loc = string( unique(Loc) );

Estimated_H = zeros(N, length(loc)); 

for i = 1:N
    a = find(pt(i) == H(:,1));
    Estimated_H(i,1) = pt(i);
    
    locs = string(Loc(a));
    
    
    j = 1;
    while j <= length(a)
        if locs(j) == loc(1)
            Estimated_H(i,2) = D( a(j));
        end 
            
        if locs(j) == loc(2)
            Estimated_H(i,3) = D( a(j));
        end 
        
        if locs(j) == loc(3)
            Estimated_H(i,4) = D( a(j));
        end 
        
        if locs(j) == loc(5)
            Estimated_H(i,5) = D( a(j));
        end 
              
        j = j + 1;
    end 
    
end 

A = [];
for i = 1:size(Estimated_H,1)
    if ismember(Estimated_H(i,1), ca_ID)
        A = [A 1];
    else 
        A = [A 0];
    end 
end 

Case = Estimated_H(find(A), :); Control= Estimated_H(find(~A), :);

fig = figure(4);
fig.Position = [15 10 1800 1500];


for i = 2:size(Case,2)
    subplot(2,3,i-1)
    [f1,xi1] = ksdensity(Case(:,i)); 
    plot(xi1,f1, 'r-','linewidth', 2);

    hold on
    [f2,xi2] = ksdensity(Control(:,i)); 
    plot(xi2,f2,'b--','linewidth', 2);
 
   % xlim([0, .0410])
    legend(["Cases", " Controls"], "fontsize", 12)
    %title(sprintf("Feature %d (p-value = %.2f) ",i, z(3,i)),'fontweight','bold','fontsize',15 )
    ylabel("Probability",'fontsize',14 )
    %xlabel(ss(i),'fontsize',14 )
    
    %a = get(gca,'XTickLabel');
    set(gca,"FontSize",15)

    grid on
    hold off
end

