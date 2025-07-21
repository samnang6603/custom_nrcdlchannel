function rhat = GetSphericalUnitVector(phi,theta)
% Get rhat
% See Equation 7.5-23 and -24
%{
For RX, the spherical unit of arrival

                [sin(theta_(n,m,ZoA))*cos(phi_(n,m,AoA))]
rhat_(rx,n,m) = |sin(theta_(n,m,ZoA))*sin(phi_(n,m,AoA))|
                [          cos(theta_(n,m,ZoA))         ]

Similarly, for departure

                [sin(theta_(n,m,ZoD))*cos(phi_(n,m,AoD))]
rhat_(tx,n,m) = |sin(theta_(n,m,ZoD))*sin(phi_(n,m,AoD))|
                [          cos(theta_(n,m,ZoD))         ]

where:
 - n is the index of the cluster
 - m is the index of the ray
 - N is the total number of cluster (L in this codebase)
 - M is the total number of ray (M in this codebase)

%}
% We can generalize this vector without creating separate ones for Tx and
% Rx. The only variables that change are theta and phi for arrival and
% departure

% permute from (L,M) (virtually (L,M,Page|Page = 1) to rearrange to
% (Page|Page = 1,L,M) to be concatenated at the end
phi = permute(phi,[3, 1, 2]);
theta = permute(theta,[3, 1, 2]);

sintheta = sind(theta); % so it's reusable, for optimization
rhat = [sintheta.*cosd(phi)
    sintheta.*sind(phi)
    cosd(theta)    ];

end