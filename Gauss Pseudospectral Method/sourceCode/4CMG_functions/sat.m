function out = sat(in, LB, UB)
% saturation function
if in > UB
    out = UB;
elseif in < LB
    out = LB;
else
    out = in;
end
end


