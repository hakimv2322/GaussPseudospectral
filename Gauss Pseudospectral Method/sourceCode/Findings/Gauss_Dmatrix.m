%--------------------------------------------------------------------------
% Gauss_Dmatrix.m
% determines approximate differentiation matrix for Gauss pseudospectral
% method with Gauss points
%--------------------------------------------------------------------------
% D = Gauss_Dmatrix(tau)
% tau: Gauss nodes, length N
%   D: differentiation matrix, dimensions Nx(N+1)
%--------------------------------------------------------------------------
% Author: Victor Hakim
%--------------------------------------------------------------------------
function D = Gauss_Dmatrix(tau)
% uses the function LagrangePoly(), function of symbolic z

N = length(tau);

if N == 9
    run('Gauss9.m')
    D = D9;
elseif N == 19
    run('Gauss19.m')
    D = D19;
elseif N == 49
    run('Gauss49.m')
    D = D49;
elseif N == 99
    run('Gauss99.m')
    D = D99;
else

    dLagr = diff(LagrangePoly(cat(1, -1, tau)));
    for it1 = 1:N
        for it2 = 0:N
            z = tau(it1);
            D(it1, it2+1) = subs(dLagr(it2+1));
        end
    end

end
end