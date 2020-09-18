%--------------------------------------------------------------------------
% Gauss_weights.m
% determines Gaussian quadrature weights using Gauss nodes
%--------------------------------------------------------------------------
% w = Gauss_weights(tau)
% tau: Gauss nodes
%   w: Gaussian quadrature weights
%--------------------------------------------------------------------------
% Author: Daniel R. Herber, Graduate Student, University of Illinois at
% Urbana-Champaign
% Date: 06/04/2015
%--------------------------------------------------------------------------
function w = Gauss_weights(tau)
% number of nodes
N = length(tau);

if N == 9
    run('Gauss9.m')
    w = w9;
elseif N == 19
    run('Gauss19.m')
    w = w19;
elseif N == 49
    run('Gauss49.m')
    w = w49;
elseif N == 99
    run('Gauss99.m')
    w = w99;
else

    % See Benson's MIT PhD thesis, page 30. Note the typo.
    [der, ~] = lepoly(N, tau);
    w = 2./((1 - tau.^2).*der.^2);

end
end