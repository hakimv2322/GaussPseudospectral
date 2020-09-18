%--------------------------------------------------------------------------
% LQR_Solution_State.m
% Calculates the optimal state for the LQR problem
% on page 126 of David Benson's PhD thesis
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
%--------------------------------------------------------------------------
function X = LQR_Solution_State(t)

    X = exp(-sqrt(2)*t) - 7.2135*(10e-7)*exp(sqrt(2)*t);
    
end