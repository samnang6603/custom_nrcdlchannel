function static_struct = GenerateStaticCDLChannel(cdl_struct,pdp_struct,ant_struct,coupling,Phi)
%GENERATESTATICCDLCHANNEL static cdl channel without Doppler effect
% TR 38.901

P = ant_struct.NumInputSignals; % Num Tx
R = ant_struct.NumOutputSignals; % Num Rx

input_size = [1 P]; % just expand dimensionality for later, that's all
M = RayModel.GetNumberOfRays;
T = 1; % 1 time-domain channel sample
D = [T 1 M P R];
L = length(info_struct.ClusterTypes);


[HLOS,rhatTxLOS,rhatRxLOS] = computeClusterGains(cdl_struct,pdp_struct,input_size,coupling,XPR,Phi,D,clusterTypes,'LOS');



end

function [Hstatic,rhat_tx,rhat_rx] = computeClusterGains(cdl_struct,pdp_struct,input_size,coupling,XPR,Phi,Dim,allClusterTypes,thisClusterType)

c0 = physconst('LightSpeed');

% get part of pdp corresponding to this cluster type
pdp = pdp_struct.PDPTable;
thisClusterIdx = strcmpi(allClusterTypes,thisClusterType);
pdpthisCluster = pdp(thisClusterIdx,:);

% get number of clusters and set up corresponding dim of Hstatic
N = size(pdpthisCluster,1);
D(2) = N;

% if no active cluster for this cluster type then return
if ~N
    return;
end

% Get column type index
[powerIdx,AoDIdx,AoAIdx,ZoDIdx,ZoAIdx] = deal(2,3,4,5,6);

% Dimension of sig in/out and number of rays
M = D(3);
S = D(4);
U = D(5);

% Initialize channel matrix and spherical unit vector rhat
Hstatic = zeros(D,'like',1i);
rhat_tx = zeros([3 D(2:3)]);
rhat_rx = zeros([3 D(2:3)]);

% Convert cluster power to linear
P = 10.^(pdpthisCluster(:,powerIdx).'/10);

% Extract per cluster param  
% Custom CDL params To be implemented later
AngleSpreads = pdp_struct.AngleSpreads;
[C_ASD,C_ASA,C_ZSD,C_ZSA] = deal(AngleSpreads);

% Get ray coupling associated with this cluster type
coupling = coupling(thisClusterIdx,:,:);

% Get XPR associated with this cluster type
xprtmp = XPR;
XPR = repmat(xprtmp,length(clusterTypes),M); % expand into XPR Map
XPR = XPR(thisClusterIdx,:,:);

% Get Phi matrix associated with this cluster type
Phi = Phi(thisClusterIdx,:,:);

% Apply ray offset angles to AoA & AoD angles 
% ray offset in 7.5 Table 7.5-3
ray_offset_alpha = kron([0.0447, 0.1413, 0.2492, 0.3715, 0.5129, 0.6797, 0.8844, 1.1481, 1.5195, 2.1551],[1, -1]); % size 1-by-M

% Only for this LOS cluster
if (strcmpi(thisClusterType,'LOS'))
    ray_offset_alpha(:) = 0; % no ray offseting for LOS path (only first ray will be used anyway)
end

% 7.5 Equation 7.5-13 for AoA and AoD
phi_AoA = pdpthisCluster(:,AoAIdx) + C_ASA*ray_offset_alpha; % [N M]
phi_AoD = pdpthisCluster(:,AoDIdx) + C_ASD*ray_offset_alpha; % [N M]

% 7.5 Equation 7.5-18 for ZoA and ZoD
theta_ZoA = pdpthisCluster(:,ZoAIdx) + C_ZSA*ray_offset_alpha;
theta_ZoD = pdpthisCluster(:,ZoDIdx) + C_ZSD*ray_offset_alpha;

% Angle Scaling section 7.7.5.1
% TO BE IMPLEMENTED LATER


% Wrap azimuth angles to [-180, 180]
% wrap modulo-360 of 180 and center at 0 by -180
phi_AoA = mod(phi_AoA + 180,360) - 180;
phi_AoD = mod(phi_AoD + 180,360) - 180;

% Wrap zenith angles to [0,360] and map [180,360] to [180,0]
theta_ZoA = wrapZenithAngles(theta_ZoA);
theta_ZoD = wrapZenithAngles(theta_ZoD);

% Rearrange ray coupling
%
%       (Sketch from MATLAB)
%           1
%        .-----> AoA
%       /
%     AoD
%       \   3          2
%        '-----> ZoD -----> ZoA
%     
% AoD to AoA coupling
phi_AOA = phi_AOA(coupling(:,:,1));
% AoD to ZoD coupling
theta_ZOD = theta_ZOD(coupling(:,:,3));
% ZoD to ZoA coupling
theta_ZOA = theta_ZOA(coupling(:,:,2));
% rearrange ZoA in AoD order
theta_ZOA = theta_ZOA(coupling(:,:,3));

% Get the spherical unit vectors of departure for each cluster & ray 
rhat_tx = getSphericalUnitVector(phi_AoD,theta_ZoD);

% Get transmit antenna/subarray location vector
if isempty(ant_struct.SubarrayPositions)
    txRadiatorPositions = ant_struct.ElementPositions;
else
    txRadiatorPositions = ant_struct.SubarrayPositions;
end

% Carrier wavelength
lambda_0 = c0/cdl_struct.CarrierFrequency;

% Calculate the location vector dbar of Tx
dbar_tx = reshape(txRadiatorPositions,3,[])*lambda_0 + ...
    ant_struct.TransmitAntennaArray.Positions;

% Allocation for Tx field term and location term
transmitFieldTerm = complex(zeros(2,N*M,S));
transmitLocationTerm = complex(zeros(N*M,S),zeros(N*M,S));


end

function theta = wrapZenithAngles(theta)
theta = mod(theta,360);
theta(theta>180) = 360 - theta(theta>180);
end

function rhat = getSphericalUnitVector(phi,theta)
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
rhat = [sintheta.*cosd(phi)...
        sintheta.*sind(phi)...
            cosd(theta)    ];

end

