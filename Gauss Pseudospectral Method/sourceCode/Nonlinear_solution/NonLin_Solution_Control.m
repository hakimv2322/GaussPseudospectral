%--------------------------------------------------------------------------
% NonLin_Solution_Control.m
% Calculates the optimal control for the nonlinear example
% on page 135 of David Benson's PhD thesis
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
%--------------------------------------------------------------------------
function X = NonLin_Solution_Control(t)

    X = -3.4122*exp(-sqrt(2)*t) + 3.5138*(10e-5)*exp(sqrt(2)*t);
    
end