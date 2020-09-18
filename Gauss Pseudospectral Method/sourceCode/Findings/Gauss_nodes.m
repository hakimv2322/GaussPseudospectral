%--------------------------------------------------------------------------
% Gauss_nodes.m
% determines Gauss nodes
%--------------------------------------------------------------------------
% tau = Gauss_nodes(N)
%   N: number of nodes minus 1, should be an integer greater than 0
% tau: Gauss nodes
%--------------------------------------------------------------------------
% Examples:
% tau = Gauss_nodes(1)
% 0
% tau = Gauss_nodes(2)
% -0.57735 0.57735
% tau = Gauss_nodes(3)
% -0.77460 0 0.77460
%--------------------------------------------------------------------------
% Contributor: Victor Hakim
%--------------------------------------------------------------------------
function tau = Gauss_nodes(N)

if N == 9
    run('Gauss9.m')
    tau = tau9;
elseif N == 19
    run('Gauss19.m')
    tau = tau19;
elseif N == 49
    run('Gauss49.m')
    tau = tau49;
elseif N == 99
    run('Gauss99.m')
    tau = tau99;
else

    syms z;
    tau = sort(double(vpasolve(legendreP(N,z) == 0)));

end
end