%--------------------------------------------------------------------------
% BangBang_Solution_States.m
% Calculates the optimal states for the Bang-Bang problem
%--------------------------------------------------------------------------
% See page 141 of Benson's PhD thesis
%--------------------------------------------------------------------------
% Primary contributor: Victor Hakim
%--------------------------------------------------------------------------
function [X1, X2] = BangBang_Solution_States(T)

    t1 = 5.34520788;
    tf = 2*t1 - 3;
    
    X1 = zeros(1,length(T));
    X2 = zeros(1,length(T));
    
    for it = 1:length(T)
        t = T(it);
        if t < t1
            X1(it) = -0.5*t^2 + 3*t + 1;
            X2(it) = -t + 3;
        else
            X1(it) = 0.5*t^2 - tf*t + 0.5*tf^2;
            X2(it) = t - tf;
        end
    end

end