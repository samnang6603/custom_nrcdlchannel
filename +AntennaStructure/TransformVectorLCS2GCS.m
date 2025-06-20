function Y = TransformVectorLCS2GCS(orientation)
% TRANSFORMVECTORLCS2GCS  Transform vector from Local to Global Coordinate 
% System (LCS → GCS)
% Based on 3GPP TR 38.901 Section 7.1.3 and TR 36.873 Section 5.1.3
%
% INPUT:
%   orientation = [bearing (deg), downtilt (deg), slant (deg)]
%
% 3GPP defines the composite rotation matrix as:
%     R = Rz(alpha) * Ry(beta) * Rx(gamma)
% where:
%   alpha = bearing (rotation about Z-axis)
%   beta  = downtilt (rotation about Y-axis)
%   gamma = slant   (rotation about X-axis)
%
% ******************************!NOTE!******************************
% 3GPP describes this R as a rotation that **maps a vector in the GCS into 
% the LCS**. In other words, they interpret R as rotating the *coordinate 
% frame*.
%
% BUT: This implementation interprets R as an **active rotation** —
% that is, rotating the *vector itself* from LCS to its equivalent in GCS.
%
% This is mathematically valid because R is an orthogonal matrix:
%     R⁻¹ = Rᵗ ⇒ inverse and transpose yield the same result.
%
% Therefore:
%   - 3GPP’s "frame rotation" (passive) using R
%   - is equivalent to this function’s "vector rotation" (active) using R
%   - just applied in the opposite conceptual sense.

c = cosd(orientation); % [cos(alpha) cos(beta) cos(gamma)]
ca = c(1); 
cb = c(2); 
cg = c(3);

s = sind(orientation); % [sin(alpha) sin(beta) sin(gamma)]
sa = s(1); 
sb = s(2); 
sg = s(3);

Y = zeros(3);
Y(1,1) = ca.*cb;
Y(1,2) = ca*sb*sg - sa*cg;
Y(1,3) = ca*sb*cg + sa*sg;
Y(2,1) = sa*cb;
Y(2,2) = sa*sb*sg + ca*cg;
Y(2,3) = sa*sb*cg - ca*sg;
Y(3,1) = -sb;
Y(3,2) = cb*sg;
Y(3,3) = cb*cg;
end