function [H,sampleTimes,cdl_struct] = Iterate(cdl_struct)
%Ierate performs one iteration of the CDL channel and update the channel 
% gain. Returns updated struct 'cdl_struct' and channel gain 'H'.

% Get parameter structures
info_struct = cdl_struct.PDPInfo;
pdp_struct = cdl_struct.PDP;
ant_struct.TransmitAntennaArray = cdl_struct.TransmitAntennaArray;
ant_struct.NumInputSignals = prod(cdl_struct.TransmitAntennaArray.Size);
ant_struct.ReceiveAntennaArray = cdl_struct.ReceiveAntennaArray;
ant_struct.NumOutputSignals = prod(cdl_struct.ReceiveAntennaArray.Size);

% Generate initial phase
[phi,cdl_struct] = ChannelPhase.GenerateInitialPhases(cdl_struct,...
    ant_struct,info_struct);

% Compute ray coupling
coupling = RayModel.ComputeRayCoupling(cdl_struct,info_struct);

% Calculate dual mobility scatterer variables
[scattererStates,scattererSpeeds] = DualMobility.ComputeDualMobilityScattererVariables( ...
    cdl_struct,info_struct);

% Generate static CDL channel
static_struct = ChannelCoefficient.GenerateStaticCDLChannel(cdl_struct,...
    pdp_struct,ant_struct,info_struct,coupling,phi);

% Generate time varying CDL channel
[H,sampleTimes] = ChannelCoefficient.GenerateTimeVaryingCDLChannel(cdl_struct,...
    ant_struct,static_struct,scattererStates,scattererSpeeds);


end