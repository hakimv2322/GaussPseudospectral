%--------------------------------------------------------------------------
% NonLin_Solution_State.m
% Calculates the optimal state for the nonlinear example
% on page 135 of David Benson's PhD thesis
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
%--------------------------------------------------------------------------
function X = NonLin_Solution_State(t)

    X = (1.4134*exp(-sqrt(2)*t) + 8.4831*(10e-5)*exp(sqrt(2)*t)).^2;
    
end