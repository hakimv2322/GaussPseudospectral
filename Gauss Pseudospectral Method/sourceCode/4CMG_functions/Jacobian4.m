function A = Jacobian4(h, beta, delta)
d1 = delta(1); d2 = delta(2);
d3 = delta(3); d4 = delta(4);
cb = cosd(beta); sb = sind(beta);
c1 = cosd(d1); s1 = sind(d1);
c2 = cosd(d2); s2 = sind(d2);
c3 = cosd(d3); s3 = sind(d3);
c4 = cosd(d4); s4 = sind(d4);
A1 = [-cb*c1, s2, cb*c3, -s4];
A2 = [-s1, -cb*c2, s3, cb*c4];
A3 = [sb*c1, sb*c2, sb*c3, sb*c4];
A = h*[A1; A2; A3];
end