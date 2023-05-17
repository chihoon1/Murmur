function [f, alfa, rst] = mfstriangle(data, L, wf, q, m1, m2, md);
% [f, alfa, rst] = mfstriangle(data, q, m1, m2, md);
%   Multifractal analysis - find two kinds of slopes which are left/right
%   tangent at f(alfa) = -0.2 and two slopes with hurst exponent term
%   Input : data -> time series, q -> range of powers, m1 and m2 -> levels for
%   slopes: wavelets m1(=6) and m2(=12), md = show or not, default = 0 no
%   show.
%   Output : f(y) of multifractal, alfa is alfa value(x), rst consists of
%   [Hurst exponent, Left(tanget), Right(tangent), Left slope with Hurst,
%   Right slope with Hurst, bandwidth, left alfa point, right alfa point at
%   f = -0.2.]
%   calls the function "mfspectra.m"

if nargin == 6, md = 0; end;  % No show unless md == 1

[f, alfa] = mfspectra(data, L, wf, q, m1, m2);   hs = alfa(1,f == 0);  
ind = [];
for i = 1 : length(hs)
    ind = [ind sum(alfa >= hs(1,i))];
end
mm1 = min(ind);     mm2 = max(ind);         % split data into 2 sets
alfa1 = alfa(1, 1:mm1);     f1 = f(1,1:mm1);
alfa2 = alfa(1, mm2 : end); f2 = f(1, mm2:end);

%   Right : find the right alfa such that f(alfa) = -0.2
pt1 = sum(f1 <= -0.2);
if f1(1, pt1) ~= -0.2
    a1 = (f1(1,pt1+1) - f1(1,pt1))/(alfa1(1, pt1+1) - alfa1(1,pt1));
    b1 = f1(1, pt1) - a1*alfa1(1,pt1);
    pt_R = (-0.2 - b1)/a1;
else
    pt_R = alfa1(1,pt1);
    a1 = (f1(1,pt1+1) - f1(1,pt1))/(alfa1(1, pt1+1) - alfa1(1,pt1));
    b1 = f1(1, pt1) - a1*alfa1(1,pt1);
end
bd = 0.2;
x1 = [pt_R-bd pt_R  pt_R+bd];
y1 = a1*x1 + b1;


%   Left : find left alfa such that f(alfa) = -0.2
pt2 = sum(f2 >= -0.2);
if f2(1, pt2) ~= -0.2
    a2 = (f2(1,pt2+1) - f2(1,pt2))/(alfa2(1, pt2+1) - alfa2(1,pt2));
    b2 = f2(1, pt2) - a2*alfa2(1,pt2);
    pt_L = (-0.2 - b2)/a2;
else
    pt_L = alfa2(1, pt2);
    a2 = (f2(1,pt2+1) - f2(1,pt2))/(alfa2(1, pt2+1) - alfa2(1,pt2));
    b2 = f2(1, pt2) - a2*alfa2(1,pt2);
end
x2 = [pt_L-bd pt_L pt_L+bd];
y2 = a2*x2 + b2;


H = mean(alfa(1,f == 0));       
RT = a1;                LT = a2;        % two tangents at f(alfa) = -0.2
RS = -0.2/(pt_R - H);   LS = -0.2/(pt_L - H);   % slopes with mode
B = pt_R - pt_L;                        % Broadness
rst = [H LS RS LT RT B pt_L pt_R];

%   Graphical representation
if md == 1
    %figure
    ld = 3;   mk = 8;     fs = 17;
    plot(alfa, f, 'linewidth', ld);
    hold on
    plot(x1, y1, 'r', 'linewidth', ld);
    text(pt_R+0.075, -0.15, num2cell(RT), 'fontsize', fs);
    
    plot(x2, y2, 'r', 'linewidth', ld);
    text(pt_L-0.22, -0.15, num2cell(LT), 'fontsize', fs);
    
    plot([pt_L pt_R], [-0.2 -0.2], 'ro', 'markersize', mk);
    plot(H, 0, 'ro', 'markersize', mk);
    text(H-0.13, 0.1, num2cell(round(abs(H),2)), 'fontsize', fs);
    
    plot([pt_L H], [-0.2 0], 'g', 'linewidth', ld);
    text(pt_L+0.51, -0.175, num2cell(round(LS,2)), 'fontsize', fs);
    
    plot([pt_R H], [-0.2 0], 'g', 'linewidth', ld);
    text(pt_R-0.5, -0.175, num2cell(round(RS,2)), 'fontsize', fs);
    
    plot([pt_L-0.2 pt_R+0.2], [-0.2 -0.2], '-.c', 'linewidth', ld);
    text(H, -0.25, num2cell(round(pt_R-pt_L,2)), 'fontsize', fs);
%     axis([pt_L-0.4 pt_R+0.4 -0.25 0]);    
    % axis([0 2 -0.2 0]);
    hold off
    xlabel('\alpha(q)', 'fontsize', 20);    ylabel('f(\alpha(q))', 'fontsize', 20);
else
end




