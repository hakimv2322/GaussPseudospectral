%--------------------------------------------------------------------------
% BangBang_Solution_Control.m
% Calculates the optimal control for the Bang-Bang problem
%--------------------------------------------------------------------------
% See page 141 of Benson's PhD thesis
%--------------------------------------------------------------------------
% Primary contributor: Victor Hakim
%--------------------------------------------------------------------------
function U = BangBang_Solution_Control(T)

    t1 = 5.34520788;
    
    U = zeros(1,length(T));
    
    for it = 1:length(T)
        t = T(it);
        if t < t1
            U(it) = -1;
        else
            U(it) = 1;
        end
    end

end