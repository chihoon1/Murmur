close all; clear all; clc
addpath('./MatlabFunctions/');

randn('seed', 1)

lw = 2.5; set(0, 'DefaultAxesFontSize', 16);fs = 15;msize = 10;

dirName = sprintf('./data/training_data');             %# folder path
files = dir( fullfile(dirName,'*.wav') );   %# list all *.xyz files
files = {files.name}';                      %'# file names
nfi = numel(files);

c1 = [3 10]; c2 = [2 9]; c3 = [1 5];  

% c1
% a1 = [1, 213, 480, 577,722, 1489, 1641,1899, 2049, 2127, 2248, 2496, 2606, 2767, 2797,2864, 3022];
% b1 = [39, 237,575, 606, 732, 1509, 1680, 1998, 2102, 2178, 2361, 2577, 2764, 2770, 2854, 2989, 3163 ];
% 
% % c2
% a2 = [40, 238, 605, 733, 1510, 1608, 1681, 1725, 2200, 2362 ];
% b2 = [212, 479, 721, 1488, 1596, 1640, 1709,1898, 2247, 2496 ];
% 
% % c3
% a3 = [576, 1597, 1710, 1998, 2103, 2179, 2578, 2771,2855, 2990];
% b3 = [576, 1607, 1724, 2048, 2126, 2119, 2605, 2796, 2863,3021];

a1 = [366, 1417, 1418, 1448, 1587, 2334, 2447];% [3,10]
b1 = [1488];% [2,9]

L = 5;  q = -5:.1:6;
filt= [-0.075765714789341  -0.029635527645954   0.497618667632458 ...
          0.803738751805216   0.297857795605542  -0.099219543576935 ...
        -0.012603967262261   0.032223100604071];

H = zeros(nfi, 9); Loc = cell(nfi, 1);
Std_H  = zeros(nfi, 2);

LL = zeros(nfi, 2);  W = [];


for i  = 1:nfi
            % extract patient ID
            newStr = split(files{i},[" ","_","."]);
            Patient_ID = str2num(newStr{1});
            Murmer_Location = newStr{2};

            % Read audio file
            fname = fullfile(dirName,files{i});
            [data, fs] = audioread(fname);

            data = (data - mean(data))/ max(abs(data));

            LL(i,:) = [Patient_ID length(data)];
            N=length(data);  LN= log2(N); NN=floor(LN);

            y= data(1:2^NN);

            H(i,1) = Patient_ID;
            
            if ismember(i,a1)
                [a, b, c] =mfstriangle(y, 1, filt, q, 3, 10, 1);
            elseif ismember(i,b1)
                [a, b, c] =mfstriangle(y, 1, filt, q, 2, 9, 1);
            else 
                [a, b, c] =mfstriangle(y, 1, filt, q, 2, NN-2, 1);
            end
            H(i,2:end) = c;


            Loc{i} =  Murmer_Location;

end 

% separate data by class labels (murmur present, absent, unknown)
D = H(:,7);
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

for i = 1:7
    D = H(:, i+1);
    X_control = D(CO); X_case = D(CA);
    
    subplot(4,2,i)
    [f1,xi1] = ksdensity(X_case); 
    plot(xi1,f1, 'r-','linewidth', 2);

    hold on
    [f2,xi2] = ksdensity(X_control); 
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
%saveas(fig,'Multifractal.png')

%% Save into .csv files
% 
%Window_Mf.Case = H(CA,:); Window_Mf.Control = H(CO,:);
%save('Window_Mf.mat','Window_Mf')

writematrix( H(CA,:),'./case_features/Window_Mf_Case.csv');
writematrix( H(CO,:),'./control_features/Window_Mf_Control.csv');
% writematrix( H(Un,:),'Window_Mf_Un.csv');

%% 

D = H(:,7);
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

