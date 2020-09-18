% Makes a family of symbolic Lagrange interpolation polynomials,
% interpolated over a 1-D array tau.
% Function of symbolic z.
% Author: Victor Hakim

function Lagr = LagrangePoly(tau);
n = length(tau);
syms z
for it1 = 1:n
    dummy = 1;
    for it2 = 1:n
        if it2 == it1
            continue
        end
        dummy = dummy*(z - tau(it2))/(tau(it1) - tau(it2));
    end
    Lagr(it1) = dummy;
end
end