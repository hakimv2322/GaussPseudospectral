% Multiply two quaternions together.

function quat = quatMult(qA, qB)
% Assumes the scalar component is q(4)
% Note that multiplication is not communitive.

quat = qA;
qAv = qA(1:3); qBv = qB(1:3);
quat(4) = qA(4)*qB(4) - dot(qAv, qBv);
quat(1:3) = qA(4)*qBv + qB(4)*qAv + cross(qAv, qBv);

end

