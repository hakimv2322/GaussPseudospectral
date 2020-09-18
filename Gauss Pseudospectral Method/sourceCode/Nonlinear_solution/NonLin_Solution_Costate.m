%--------------------------------------------------------------------------
% NonLin_Solution_Costate.m
% Calculates the costate for the nonlinear example
% on page 135 of David Benson's PhD thesis
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
%--------------------------------------------------------------------------
function X = NonLin_Solution_Costate(t)

    X = -NonLin_Solution_Control(t)./(2*sqrt(NonLin_Solution_State(t)));
    
end