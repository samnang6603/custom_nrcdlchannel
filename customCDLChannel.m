function cdl = customCDLChannel(profile)
if nargin == 0
    profile = 'CDL-A';
end
%% Create CDL structure
cdl = struct();

%% Predefined Delay Profile
cdl.DelayProfile = profile;
cdl.DelaySpread = 3e-8;
cdl.KFactorScaling = false;
cdl.KFactor = 9.0;
cdl.AngleScaling = false;
cdl.MeanAngles = [0.0 0.0 0.0 0.0];
%% Custom Delay Profile

% cdl.PathDelays = 0;
% cdl.AveragePathGains = 0;
% cdl.AnglesAoD = 0;
% cdl.AnglesAoA = 0;
% cdl.AnglesZoD = 0;
% cdl.AnglesZoA = 0;
% cdl.HasLOSCluster = false;
% cdl.AngleSpreads = [5 11 3 3];
% cdl.MeanAngles = [0, 0, 0, 0];
%%
cdl.RayCoupling = 'Random';
cdl.InitialPhases = 'Random';
%% Antenna Array

antTx_struct.Size = [2 2 2 1 1];
antTx_struct.ElementSpacing = [0.5 0.5 1.0 1.0];
antTx_struct.PolarizationAngles = [45 -45];
antTx_struct.Orientation = [0; 0; 0];
antTx_struct.Element = '38.901';
antTx_struct.PolarizationModel = 'Model-2';

antRx_struct.Size = [1 1 2 1 1];
antRx_struct.ElementSpacing = [0.5 0.5 0.5 0.5];
antRx_struct.PolarizationAngles = [0 90];
antRx_struct.Orientation = [0; 0; 0];
antRx_struct.Element = 'isotropic';
antRx_struct.PolarizationModel = 'Model-2';

cdl.TransmitAntennaArray = antTx_struct;
cdl.TransmitArrayOrientation = [0; 0; 0];
cdl.ReceiveAntennaArray = antRx_struct;
cdl.ReceiveArrayOrientation = [0; 0; 0];
%cdl = hArrayGeometry(cdl,NTx,NRx);
cdl.CarrierFrequency = 4e+9; % Start with 4Ghz

%% Mobility
cdl.MaximumDopplerShift = 5;
cdl.UTDirectionOfTravel = [0; 90]; % User terminal direction of travel
cdl.MovingScattererProportion= 0.2;
cdl.MaximumScattererSpeed = 5;

%% Channel Control
cdl.SampleRate = 30720000;
cdl.InitialTime = 0.0; % Time offset of fading process in seconds
cdl.SampleDensity = 64; % Number off time samples per half wavelength
cdl.RandomStream = 'mt19937ar with seed';
cdl.Seed = 73;
cdl.NormalizeChannelOutputs = true;
cdl.NormalizePathGains = true;
cdl.ChannelResponseOutput = 'path-gains';
cdl.ChannelFiltering = false;
cdl.NumTimesSamples = 30720;
cdl.OutputDataType = 'double';

%% Channel Filtering Setting
cdl.ChannelFilter.SampleRate = cdl.SampleRate;
cdl.ChannelFilter.ChannelFilterDelay = 7;
cdl.ChannelFilter.StopbandAttenuation = 70; % dB
cdl.ChannelFilter.MaxFractionalDelayError = 0.01;

end