% Find the inverse of a quaternion.

function quat = quatInv(q)
% Assumes the scalar component is q(4)

qMagSq = q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2;
quat = q/qMagSq;
quat(1:3) = -quat(1:3);

end

