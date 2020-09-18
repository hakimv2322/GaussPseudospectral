function hTotal = angVector3(h, beta, delta)
d1 = delta(1); d2 = delta(2); d3 = delta(3);
cb = cosd(beta); sb = sind(beta);
c1 = cosd(d1); s1 = sind(d1);
c2 = cosd(d2); s2 = sind(d2);
c3 = cosd(d3); s3 = sind(d3);
h1 = [-cb*s1; c1; sb*s1];
h2 = [-c2; -cb*s2; sb*s2];
h3 = [cb*s3; -c3; sb*s3];
hTotal = h*(h1 + h2 + h3);
end