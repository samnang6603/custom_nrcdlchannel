function [H,sampleTimes,cdl_struct] = GenerateCDLChannel(cdl_struct)
% Generates 5G NR Clustered Delay Line (CDL) Channel. CDL simulates how 
% radio signals propagate, particularly when the received signal is 
% composed of multiple delayed "clusters" of multipath components. 
% Each cluster has a distinct delay, but within a cluster, the many 
% individual paths have slightly different angles of departure and arrival, 
% representing reflections and scattering from different objects in the 
% environment.

% Input:
% - cdl_struct: Structure containing CDL Paramaters needed to construct the
% channel
%
% Output:
% - H: Channel impulse response
% - sampleTimes: Sample times for channel filtering operatioons
% - cdl_struct: includes new parameters based on unique channel
% characteristics and computations

% Inspired by MATLAB nrCDLChannel()

%% Step 1: configure antenna structure
[cdl_struct,ant_struct] = AntennaStructure.ConfigureStructure(cdl_struct);

%% Step 2: configure delay profile
pdp_struct = PowerDelayProfile.CreatePowerDelayProfileStructure(cdl_struct);
cdl_struct.PDP = pdp_struct;

%% Step 3: split LOS an NLOS cluster and sub-clustering
[pdp_struct,info_struct] = PowerDelayProfile.GetClusterInfo(pdp_struct);
cdl_struct.PDPInfo = info_struct;

%% Step 4: generate initial phase
cdl_struct = initializeRandomStream(cdl_struct);
[phi,cdl_struct] = ChannelPhase.GenerateInitialPhases( ...
    cdl_struct,ant_struct,info_struct);

%% Step 5: compute ray coupling
coupling = RayModel.ComputeRayCoupling(cdl_struct,info_struct);

%% Step 6: calculate dual mobility scatterer variables
[scattererStates,scattererSpeeds] = DualMobility.ComputeDualMobilityScattererVariables( ...
    cdl_struct,info_struct);
cdl_struct.Scatterer.States = scattererStates;
cdl_struct.Scatterer.Speeds = scattererSpeeds;

%% Step 7: generate static CDL Channel
static_struct = ChannelCoefficient.GenerateStaticCDLChannel( ...
    cdl_struct,pdp_struct,ant_struct,info_struct,coupling,phi);
cdl_struct.StaticChannel = static_struct;

%% Step 8: generate time varying CDL Channel
[H,sampleTimes] = ChannelCoefficient.GenerateTimeVaryingCDLChannel(cdl_struct,...
    ant_struct,static_struct,scattererStates,scattererSpeeds);
cdl_struct = ChannelProcess.AdvanceTime(cdl_struct);

end

%% Local helper fcn
function cdl_struct = initializeRandomStream(cdl_struct)
switch lower(cdl_struct.RandomStream)
    case 'mt19937ar with seed'
        randomStream = RandStream('mt19937ar',Seed=cdl_struct.Seed);
        cdl_struct.RandomStreamObj = randomStream;
    otherwise
        rng(cdl.Seed)
        cdl_struct.RandomStreamObj = [];
end
end

