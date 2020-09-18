%--------------------------------------------------------------------------
% LQR_Solution_Costate.m
% Calculates the costate for the LQR problem
% on page 126 of David Benson's PhD thesis
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
%--------------------------------------------------------------------------
function X = LQR_Solution_Costate(t)

    X = 2.4124*exp(-sqrt(2)*t) - 2.9879*(10e-7)*exp(sqrt(2)*t);
    
end