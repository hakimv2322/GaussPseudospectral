%--------------------------------------------------------------------------
% BangBang_Solution_Costates.m
% Calculates the costates for the Bang-Bang problem
%--------------------------------------------------------------------------
% See page 141 of Benson's PhD thesis
%--------------------------------------------------------------------------
% Primary contributor: Victor Hakim
%--------------------------------------------------------------------------
function [X1, X2] = BangBang_Solution_Costates(T)

    t1 = 5.34520788;
    tf = 2*t1 - 3;
    
    c1 = 1/(tf - t1);
    c2 = c1*t1;
    
    X1 = c1*ones(size(T));
    X2 = -c1*T + c2;

end