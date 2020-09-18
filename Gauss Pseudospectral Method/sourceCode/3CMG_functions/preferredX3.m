function delta = preferredX3(h, beta)
dt = 1/1000;
T = [1; 0; 0];
d1 = -88; % nominal -90
d2 = 180; % nominal 180
d3 = 88;  % nominal 90
delta = [d1; d2; d3];
while norm(angVector3(h, beta, delta)) > 0.00001*h
    Ddelta = inv(Jacobian3(h, beta, delta))*T;
    Ddelta = Ddelta*180/pi;
    delta = delta - dt*Ddelta;
%     norm(angVector3(h, beta, [d1; d2; d3]))
end
end
% Result:
% delta = [56.4; 180; -56.4];
% This corresponds to delta3 from Table 7 of Lee et al. (2018)

