function [H,sampleTimes,cdl_struct] = GenerateCDLChannel(cdl_struct)
% Generate NR CDL Channel, 
% inspired by MATLAB nrCDLChannel()

% Step 1: configure antenna structure
ant_struct = AntennaStructure.ConfigureStructure(cdl_struct);

% Step 2: configure delay profile
pdp_struct = PowerDelayProfile.CreatePowerDelayProfileStructure(cdl_struct);
cdl_struct.PDP = pdp_struct;

% Step 3: split LOS an NLOS cluster and sub-clustering
[pdp_struct,info_struct] = PowerDelayProfile.GetClusterInfo(pdp_struct);
cdl_struct.PDPInfo = info_struct;

% Step 4: generate initial phase
[Phi,cdl_struct] = ChannelPhase.GenerateInitialPhases( ...
    cdl_struct,ant_struct,info_struct);

% Step 5: compute ray coupling
coupling = RayModel.ComputeRayCoupling(cdl_struct,info_struct);

% Step 6: calcualte dual mobility scatterer variables
[scattererStates,scattererSpeeds] = DualMobility.ComputeDualMobilityScattererVariables( ...
    cdl_struct,info_struct);

% Step 7: generate static CDL Channel
static_struct = ChannelCoefficient.GenerateStaticCDLChannel( ...
    cdl_struct,pdp_struct,ant_struct,info_struct,coupling,Phi);

% Step 8: generate time varying CDL Channel
[H,sampleTimes] = ChannelCoefficient.GenerateTimeVaryingCDLChannel(cdl_struct,...
    ant_struct,static_struct,scattererStates,scattererSpeeds);

end