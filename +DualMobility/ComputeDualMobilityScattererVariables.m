function [scattererStates,scattererSpeeds] = ComputeDualMobilityScattererVariables(cdl_struct,info_struct)
%COMPUTEDUALMOBILITYSCATTERERVARIABLES Get moving scatterer variables

clusterTypes = info_struct.ClusterTypes;

% Get the LOS cluster
LOS = strcmpi(clusterTypes,'LOS'); 

% Get the NLOS subcluster indices
splitNLOS = find(strcmpi(clusterTypes,'SubclusteredNLOS'));

% Scatterer proportion dictates how many % of moving scatterers relative to
% all clusters
p = cdl_struct.MovingScattererProportion;

% Scatterers are FIXED in CDL. 
% Tx/Rx motion creates time-varying effects: Doppler, delay, angle 
% evolution. The term "moving scatterer" = lazy abstraction of relative 
% motion from Rx’s PoV.
% You model Doppler, delay slope, etc. as if the rays came from 
% moving points.
% But no one is updating scatterer positions frame-by-frame.
% So here, v actually represents the moving Tx speed that causes direction
% shift when arrived at Rx
v_scatt = cdl_struct.MaximumScattererSpeed;

if isscalar(cdl_struct.MaximumDopplerShift)
    v_scatt = 0;
end

M = RayModel.GetNumberOfRays;
N = length(clusterTypes);

randomStream = cdl_struct.RandomStreamObj;

if isempty(randomStream) % global stream
    randAlpha = rand(M,N);
    randD = rand(M,N);
else
    randAlpha = rand(randomStream,M,N);
    randD = rand(randomStream,M,N);
end

% Find out how many moving scatterers
scattererStates = double(randAlpha < p);

% For LOS cluster, there is no moving scatterers
scattererStates(LOS,:) = 0;

% Generate RV D (size [N M]) between -v_scatt to v_scatt
% randD is (0,1) so (0,1)*2-1 = (-1,1). (-1,1)*v_scatt = (-v_scatt,v_scatt)
scattererSpeeds = (randD*2 - 1)*v_scatt;

% To be implemented later, available subcluster scenario

