% Find the conjugate of a quaternion.

function quat = quatConj(q)
% Assumes the scalar component is q(4)

quat = q;
quat(1:3) = -quat(1:3);

end

