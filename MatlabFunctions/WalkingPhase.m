function [a_st, a_sw, VGRF]  = WalkingPhase(t, A, F)
% t - time 
% A - verticle ground reaction force (VGRF)
% F - VGRF valude decides stance and swing phase 

a_st = []; % Swing phase 
a_sw = []; % Stance Phase
VGRF = []; % 

j = 1;
while j < length(A)
    
    % stance phase 
    b = [];
    while j <= length(A)  &&  A(j) > F
        b  = [b j];
        j = j + 1;     
    end 
    if length(b) >=2
        a_st = [a_st t(b(end)) - t( b(1) )];
    end 
    
    if length(b) >= 10
        GRF = A(b);
        VGRF = [VGRF; [ max(GRF(1:5)) max(GRF(end-5:end))]];
    end 
    
    if j == length(A) + 1
        break
    end 
   
    % swing phase 
   c = [];
    while j <= length(A) && A(j) <= F 
        c  = [c j];
        j = j + 1;
    end 
    if length(c) >= 2
        a_sw = [a_sw t( c(end) ) - t(c(1) )];
    end 
end 
