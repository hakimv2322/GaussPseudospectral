function output = Limiter(A, err, J, slewMax, deltaMax, k, c)
% deltaMax and slewMax must be in rad/s
deltaTemp = [0; 0; 0; 0];
maxArray = [0; 0; 0];
output = [0; 0; 0];

for it1 = 1:2
    for it2 = 1:2
        for it3 = 1:2
            for it4 = 1:2
                deltaTemp(1) = deltaMax*(-1)^(it4-1);
                deltaTemp(2) = deltaMax*(-1)^(it3-1);
                deltaTemp(3) = deltaMax*(-1)^(it2-1);
                deltaTemp(4) = deltaMax*(-1)^(it1-1);
                htemp = A*deltaTemp;
                for it5 = 1:3
                    aTemp = abs(htemp(it5)/J(it5,it5));
                    maxArray(it5) = max(maxArray(it5),aTemp);
                end
            end
        end
    end
end

for it = 1:3
%     output(it) = 0.5*(c/k)*min(sqrt(4*0.4*maxArray(it)*abs(err(it))),slewMax);
    output(it) = 0.5*(c/k)*min(sqrt(4*0.4*(1000/J(it,it))*abs(err(it))),slewMax); % temporary
end
end


