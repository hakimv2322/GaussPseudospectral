function hTotal = angVector4(h, beta, delta)
d1 = delta(1); d2 = delta(2);
d3 = delta(3); d4 = delta(4);
cb = cosd(beta); sb = sind(beta);
c1 = cosd(d1); s1 = sind(d1);
c2 = cosd(d2); s2 = sind(d2);
c3 = cosd(d3); s3 = sind(d3);
c4 = cosd(d4); s4 = sind(d4);
h1 = [-cb*s1; c1; sb*s1];
h2 = [-c2; -cb*s2; sb*s2];
h3 = [cb*s3; -c3; sb*s3];
h4 = [c4; cb*s4; sb*s4];
hTotal = h*(h1 + h2 + h3 + h4);
end