function [H,sampleTimes,cdl_struct] = AdvanceIteration(cdl_struct)
%AdvanceIteration performs one iteration of the CDL channel and update the channel 
% gain. Returns updated struct 'cdl_struct' and channel gain 'H'.

% Get parameter structures
info_struct = cdl_struct.PDPInfo;
pdp_struct = cdl_struct.PDP;
ant_struct.TransmitAntennaArray = cdl_struct.TransmitAntennaArray;
ant_struct.NumInputSignals = prod(cdl_struct.TransmitAntennaArray.Size);
ant_struct.ReceiveAntennaArray = cdl_struct.ReceiveAntennaArray;
ant_struct.NumOutputSignals = prod(cdl_struct.ReceiveAntennaArray.Size);
static_struct = cdl_struct.StaticChannel;
scattererStates = cdl_struct.Scatterer.States;
scattererSpeeds = cdl_struct.Scatterer.Speeds;

% Generate time varying CDL channel
[H,sampleTimes] = ChannelCoefficient.GenerateTimeVaryingCDLChannel(cdl_struct,...
    ant_struct,static_struct,scattererStates,scattererSpeeds);
cdl_struct = ChannelProcess.AdvanceTime(cdl_struct);


end