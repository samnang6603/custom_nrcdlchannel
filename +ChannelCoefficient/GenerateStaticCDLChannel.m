function static_struct = GenerateStaticCDLChannel(cdl_struct,pdp_struct,ant_struct,info_struct,coupling,Phi)
%GENERATESTATICCDLCHANNEL static cdl channel without Doppler effect
% TR 38.901

P = ant_struct.NumInputSignals; % Num Tx
R = ant_struct.NumOutputSignals; % Num Rx

clusterTypes = info_struct.ClusterTypes;

%input_size = [1 P]; % just expand dimensionality for later, that's all
M = RayModel.GetNumberOfRays;
T = 1; % 1 time-domain channel sample
D = [T 1 M P R];
L = length(clusterTypes);

% XPR
xprtmp = pdp_struct.XPR;
XPR = repmat(xprtmp,L,M); % expand into XPR Map


[HLOS,rhatTxLOS,rhatRxLOS] = computeClusterGains( ...
    cdl_struct,pdp_struct,ant_struct,...
    coupling,XPR,Phi,D,clusterTypes,'LOS');

% Subcluster to be implemented later
%[HNLOSubcluster,rhatTxSubcluster,rhatRxSubcluster]

[HNLOS,rhatTxNLOS,rhatRxNLOS] = computeClusterGains( ...
    cdl_struct,pdp_struct,ant_struct,...
    coupling,XPR,Phi,D,clusterTypes,'NLOS');

% Gather all channel matrices and unit vectors
LOSIdx = strcmpi(clusterTypes,'LOS');
% Subcluster to be implemented later
NLOSIdx = strcmpi(clusterTypes,'NLOS');

% Static channel matrices
pathGains = zeros(1,L,M,P,R,'like',1i);
pathGains(:,LOSIdx,:,:,:) = HLOS;
pathGains(:,NLOSIdx,:,:,:) = HNLOS;

% Departure ray's spherical unit vector
txSphericalUnitVectors = zeros(3,L,M);
txSphericalUnitVectors(:,LOSIdx,:) = rhatTxLOS;
txSphericalUnitVectors(:,NLOSIdx,:) = rhatTxNLOS;

% Arrival ray's spherical unit vector
rxSphericalUnitVectors = zeros(3,L,M);
rxSphericalUnitVectors(:,LOSIdx,:) = rhatRxLOS;
rxSphericalUnitVectors(:,NLOSIdx,:) = rhatRxNLOS;

%  Normalizations
if cdl_struct.NormalizeChannelOutputs
    pathGains = pathGains/sqrt(size(pathGains,5));
end

if cdl_struct.NormalizePathGains
    powIdx = 2;
    % cluster powers in TR 38.901, AKA average path gains
    avgPathGains = sum(10.^(pdp_struct.PDPTable(:,powIdx)/10)); 
    pathGains = pathGains/sqrt(avgPathGains);
end

static_struct = struct();
static_struct.StaticPathGains = pathGains;
static_struct.txSphericalUnitVectors = txSphericalUnitVectors;
static_struct.rxSphericalUnitVectors = rxSphericalUnitVectors;

end

function [Hstatic,rhat_tx,rhat_rx] = computeClusterGains( ...
    cdl_struct,pdp_struct,ant_struct,coupling,XPR,Phi,D,...
    allClusterTypes,thisClusterType)

c0 = 299792458;

% get part of pdp corresponding to this cluster type
pdp = pdp_struct.PDPTable;
thisClusterIdx = strcmpi(allClusterTypes,thisClusterType);
pdpthisCluster = pdp(thisClusterIdx,:);

% get number of this clusters and set up corresponding dim of Hstatic
N = size(pdpthisCluster,1);
D(2) = N;

% if no active cluster for this cluster type then return
if ~N
    return;
end

% Get column type index
[powerIdx,AoDIdx,AoAIdx,ZoDIdx,ZoAIdx] = deal(2,3,4,5,6);

% Dimension of sig in/out and number of rays
% Naming var M, S, U to be consistent with TR 38.901
M = D(3); 
S = D(4); % Total number of Tx
U = D(5); % Total number of Rx

% Initialize channel matrix and spherical unit vector rhat
Hstatic = zeros(D,'like',1i);
% rhat_tx = zeros([3 D(2:3)]);
% rhat_rx = zeros([3 D(2:3)]);

% Convert cluster power to linear
P = 10.^(pdpthisCluster(:,powerIdx).'/10);

% Extract per cluster param  
% Custom CDL params To be implemented later
AngleSpreads = pdp_struct.AngleSpreads;
C_ASD = AngleSpreads(1);
C_ASA = AngleSpreads(2);
C_ZSD = AngleSpreads(3);
C_ZSA = AngleSpreads(4);

% Get ray coupling associated with this cluster type
coupling = coupling(thisClusterIdx,:,:);

% Get XPR associated with this cluster type
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
phi_AoA = phi_AoA(coupling(:,:,1));
% AoD to ZoD coupling
theta_ZoD = theta_ZoD(coupling(:,:,3));
% ZoD to ZoA coupling
theta_ZoA = theta_ZoA(coupling(:,:,2));
% rearrange ZoA in AoD order
theta_ZoA = theta_ZoA(coupling(:,:,3));

% Carrier wavelength
lambda_0 = c0/cdl_struct.CarrierFrequency;

% Processing for Tx ------------------------------------------------------
% Get Tx array struct
txArray_struct = ant_struct.TransmitAntennaArray;

% Get the spherical unit vectors of departure for each cluster & ray 
rhat_tx = getSphericalUnitVector(phi_AoD,theta_ZoD);

% Get transmit antenna/subarray location vector
if isempty(txArray_struct.SubarrayPositions)
    txRadiatorPositions = txArray_struct.ElementPositions;
else
    txRadiatorPositions = txArray_struct.SubarrayPositions;
end

% Calculate the location vector dbar of Tx
dbar_tx = reshape(txRadiatorPositions,3,[])*lambda_0 + ...
    txArray_struct.Position;

% Allocation for Tx field term and location term
txFieldTerm = complex(zeros(2,N*M,S));
txLocationTerm = complex(zeros(N*M,S),zeros(N*M,S));
for s = 1:S
    % Calculate individual antenna element field
    txFieldTerm(:,:,s) = AntennaStructure.ElementPattern.ComputeFieldTerm(...
        txArray_struct,theta_ZoD,phi_AoD,s);

    % Get individual antenna element location
    txLocationTerm(:,s) = getLocationTerm(rhat_tx,dbar_tx,lambda_0,s);

end
% End Processing for Tx ---------------------------------------------------

% Processing for Rx ------------------------------------------------------
% Get Rx array struct
rxArray_struct = ant_struct.ReceiveAntennaArray;

% Get the spherical unit vectors of departure for each cluster & ray 
rhat_rx = getSphericalUnitVector(phi_AoA,theta_ZoA);

% Get transmit antenna/subarray location vector
if isempty(rxArray_struct.SubarrayPositions)
    rxRadiatorPositions = rxArray_struct.ElementPositions;
else
    rxRadiatorPositions = rxArray_struct.SubarrayPositions;
end

% Calculate the location vector dbar of Tx
dbar_rx = reshape(rxRadiatorPositions,3,[])*lambda_0 + ...
    rxArray_struct.Position;

% Allocation for Tx field term and location term
rxFieldTerm = complex(zeros(2,N*M,U));
rxLocationTerm = complex(zeros(N*M,U),zeros(N*M,U));
for u = 1:U
    % Calculate individual antenna element field
    rxFieldTerm(:,:,u) = AntennaStructure.ElementPattern.ComputeFieldTerm(...
        rxArray_struct,theta_ZoA,phi_AoA,u);

    % Get individual antenna element location
    rxLocationTerm(:,u) = getLocationTerm(rhat_rx,dbar_rx,lambda_0,u);
end
% End Processing for Tx ---------------------------------------------------

% Calculate polarization term
polTerm = calculatePolarizationMatrix(Phi,XPR,thisClusterType);

% Get active ray(s) in this clusterr
R = findActiveRays(Phi);

% for each Tx
for s = 1:S

    % for each Rx
    for u = 1:U

        % Calculate the MONSTROUS Equation 7.5-28 excluding Doppler (the
        % last term). Doppler is to be calculated in
        % GenerateTimeVaryingCDLChannel()
        % First, compute all the field terms combined. In 7.5-28, all the 
        % field terms are combined as below:
        %{ 
            (thth): theta-theta, (thph): theta-phi
            ksqrinv: 1/sqrt(kappa)

fieldTerms =
         T
[F_rx_th]   [exp(jPhi_thth)          ksqrinv*exp(jPhi_thph]   [F_tx_th]
|       | * |                                             | * |       |
[F_rx_ph]   [ksqrinv*exp(jPhi_phth)         exp(jPhi_phph)]   [F_rx_ph]
             
        %}
        % This implmentation doesn't
        % include matrix multiplication as in 7.5-28 but an expanded
        % vectorized version for optimization.
        fieldTerms = (rxFieldTerm(1,:,u).*polTerm(1,:) + rxFieldTerm(2,:,u).*polTerm(2,:)).*txFieldTerm(1,:,s) +...
                     (rxFieldTerm(1,:,u).*polTerm(3,:) + rxFieldTerm(2,:,u).*polTerm(4,:)).*txFieldTerm(2,:,s);
        fieldTerms = fieldTerms(:); % collapse

        % The, compute all the location terms combined. In 7.5-28, all the 
        % location terms are combined as below (excluding Doppler term)
        %{
                           T                         T
                    /     r_rx*dbar_rx \      /     r_tx*dbar_tx \    
                 exp|j2pi--------------| * exp|j2pi--------------|
                    \       lambda_0   /      \       lambda_0   / 
    
        %}
        locationTerms = rxLocationTerm(:,u).*txLocationTerm(:,s);

        % Now multiply fieldTerms and locationTerms
        allTerms = fieldTerms.*locationTerms;

        % Force zero to indices that correspond to unused rays for the
        % current set of clusters
        allTerms(R(:)==0) = 0;

        % Hstatic
        Hstatic(1,:,:,s,u) = reshape(allTerms,N,M);

    end

end

% Number of active rays for each cluster, |R_i|
modR_i = sum(R,2).';

% The scaling factor here is the sqrt(P/M) where M is the number of active
% ray per cluster
scaling = sqrt(P./modR_i);

% Now scale Hstatic
Hstatic = scaling.*Hstatic;

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
rhat = [sintheta.*cosd(phi)
    sintheta.*sind(phi)
    cosd(theta)    ];

end

function locTerm = getLocationTerm(rhat,dbar,lambda_0,antIdx)
% reshape rhat to combine cluster and ray dimensions into a single row
rhat = reshape(rhat,[3 numel(rhat)/3]);
locTerm = exp(1i*2*pi*rhat.'*dbar(:,antIdx)/lambda_0);
end

function polmat = calculatePolarizationMatrix(Phi,XPR,clusterType)
switch clusterType
    case 'LOS'
        % Easy one line, LOS is not complicated
        polmat = [1; 0; 0; -1].*exp(1i*Phi(:,:,1));
    case 'NLOS'
        % Reshape Phi to combine cluster and ray dims into a single
        % row
        % First bring the polarization dim to first dim 3D-"Transpose"
        tmp = permute(Phi,[3 1 2]);

        % Then reshape to concatenate everything in 1st dim
        Phi = reshape(tmp,[4 numel(Phi)/4]);

        % Compute kappa: the cross-polarization power ratio (XPR), given in dB.
        % XPR is defined as the ratio of power in co-polarization (theta-theta or phi-phi)
        % to cross-polarization (theta-phi or phi-theta). See TR 38.901 Eq 7.5-21.
        kappa = 10.^(XPR/10);

        % Reshape and permute so that kappa is a 1×N vector (flattened ray dimension),
        % matching how Phi is arranged (basically [numCluster×numRays]=N).
        % This aligns dimensions for vectorized element-wise multiplication later.
        tmp = permute(kappa,[3 1 2]);
        kappa = reshape(tmp,[1 numel(kappa)]);

        % Calculate complex exponentials of initial phases.
        % Phi is expected to be 4×N (for theta-theta, theta-phi, phi-theta, phi-phi).
        expPhi = exp(1i*Phi);  % This is already in the j-omega jungle.

        % For polarization matrix, we apply per-ray scaling based on kappa.
        % See Eq 7.5-22 — except we now vectorize the matrix form into columns for faster sim.

        % Define scalar scaling factors for each matrix element.
        % tmp12 and tmp21 are the inverse sqrt(kappa) terms, used for cross-polarized components.
        tmp11 = ones(size(kappa));           % theta-theta and phi-phi use 1
        tmp12 = 1./sqrt(kappa);            % theta-phi (tmp21) and phi-theta (tmp12) are scaled down
        tmp21 = tmp12;                       % symmetric cross-polar terms
        tmp22 = tmp11;                       % again, co-polar

        % Construct the polarization matrix terms, one column per ray:
        % Order is critical! MATLAB reshapes Phi in row-major order, done
        % in the GenerateInitialPhase() function, while 3GPP's definition 
        % is effectively column-major (see 7.5-22). So:
        %   Phi(1,:) = theta-theta
        %   Phi(2,:) = phi-theta  <-- swapped from intuition
        %   Phi(3,:) = theta-phi  <-- this is also out-of-order
        %   Phi(4,:) = phi-phi
        %
        % This is what happens when standards and MATLAB memory layouts 
        % have a turf war.
        %
        P1 = tmp11.*expPhi(1,:);   % theta-theta
        P2 = tmp21.*expPhi(3,:);   % theta-phi — careful! index 3 due to weird reshape ordering
        P3 = tmp12.*expPhi(2,:);   % phi-theta — index 2 despite sounding like 3
        P4 = tmp22.*expPhi(4,:);   % phi-phi

        % Stack columns to form a 4×N polarization matrix
        % (where each column is a ray's 2x2 matrix in column-major order):
        % [P1; P2; P3; P4] means:
        % [ theta-theta ;
        %   theta-phi ;
        %   phi-theta ;
        %   phi-phi ] per ray.
        polmat = [P1; P2; P3; P4];
end

end

function R = findActiveRays(Phi)
[N,M,~] = size(Phi);
R = ones(N,M);
R(Phi(:,:,1)==-Inf) = 0;
end