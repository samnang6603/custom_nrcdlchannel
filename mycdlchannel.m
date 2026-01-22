%% My custom CDL Channel Model, Default param same as nrCDLChannel()
clear, clc
%% Create CDL structure
cdl = struct();

%% Define NTx, NRx
NTx = 2;
NRx = 2;

%% Predefined Delay Profile
cdl.DelayProfile = 'CDL-D';
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
cdl.InitialTime = 0.0; % Initial time offset of fading process in seconds
cdl.CurrentTime = cdl.InitialTime; % Current time offset of fading process in seconds
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

%% Process
%model = makecdlchannelstructure(cdl);
tic
[H11,sampleTimes,cdl] = GenerateCDLChannel(cdl);
toc
[H12,sampleTimes,cdl] = ChannelProcess.AdvanceIteration(cdl); 

% Test with input signal filtering
% MATLAB CDL channel defaults follow 3GPP 5G NR/LTE numerology:
% - SampleRate = 30.72 MHz (standard for 15 kHz SCS in FR1/LTE)
% - NumTimeSamples = 30720 -> corresponds to 1 ms waveform
%   since 30720 / 30.72e6 = 1 ms exactly
%
% This is aligned with 1 subframe (1 ms) duration in LTE
% and maps well to OFDM-based waveform simulations.
%
% A subframe in 5G/LTE includes:
% - Multiple OFDM symbols (14 per slot)
% - All RBs spanning the carrier bandwidth
% - Structured according to the subcarrier spacing (SCS)
%
% This makes CDL response suitable for subframe-level analysis
% in PHY-layer simulations (BLER, EVM, etc.).
rng(10)
sigin = rand(30720,8,like=1i); % a 1ms by 5G standards subframe
tic
sigout = ChannelFilter.FilterCDLChannel(sigin,H,sampleTimes,cdl);
toc




%% Local Fcn Array Geometry
function cdl = hArrayGeometry(cdl,NTxAnts,NRxAnts,varargin)

    if (nargin==3)
        linkDirection = 'downlink';
    else
        linkDirection = varargin{1};
    end

    if (strcmpi(linkDirection,'downlink'))
        txArray = bsArrayGeometry(cdl.TransmitAntennaArray,NTxAnts);
        rxArray = ueArrayGeometry(cdl.ReceiveAntennaArray,NRxAnts);
    else % uplink
        txArray = ueArrayGeometry(cdl.TransmitAntennaArray,NTxAnts);
        rxArray = bsArrayGeometry(cdl.ReceiveAntennaArray,NRxAnts);
    end

    % Update CDL channel arrays configuration
    cdl.TransmitAntennaArray = txArray;
    cdl.ReceiveAntennaArray = rxArray;

    %warnIfArraySizeChanged(cdl,NTxAnts,NRxAnts,linkDirection)

end

function Array = bsArrayGeometry(Array,nBsAnts)

    % Setup the base station antenna geometry
    % Table of antenna panel array configurations
    % M=  no. of rows in each antenna panel
    % N=  no. of columns in each antenna panel
    % P=  no. of polarizations (1 or 2)
    % Mg= no. of rows in the array of panels
    % Ng= no. of columns in the array of panels
    % Row format= [M  N   P   Mg  Ng]
    antArraySizes = ...
       [1   1   1   1   1;   % 1 ants
        1   1   2   1   1;   % 2 ants
        2   1   2   1   1;   % 4 ants
        2   2   2   1   1;   % 8 ants
        2   4   2   1   1;   % 16 ants
        4   4   2   1   1;   % 32 ants
        4   4   2   1   2;   % 64 ants
        4   8   2   1   2;   % 128 ants
        4   8   2   2   2;   % 256 ants
        8   8   2   2   2;   % 512 ants
        8  16   2   2   2];  % 1024 ants
    antselected = min(1+ceil(log2(nBsAnts)),size(antArraySizes,1));
    Array.Size = antArraySizes(antselected,:);

    % Adjust element spacing to avoid panel overlaps
    Array.ElementSpacing(3) = Array.Size(1)*Array.ElementSpacing(1);
    Array.ElementSpacing(4) = Array.Size(2)*Array.ElementSpacing(2);

end

function array = ueArrayGeometry(array,nUeAnts)

    % Setup the UE antenna geometry
    if nUeAnts == 1
        % In the following settings, the number of rows in antenna array, 
        % columns in antenna array, polarizations, row array panels and the
        % columns array panels are all 1
        arraySize = ones(1,5);
    else
        % In the following settings, the no. of rows in antenna array is
        % nUeAnts/2, the no. of columns in antenna array is 1, the no.
        % of polarizations is 2, the no. of row array panels is 1 and the
        % no. of column array panels is 1. The values can be changed to
        % create alternative antenna setups
        arraySize = [ceil(nUeAnts/2),1,2,1,1];
    end
    array.Size = arraySize;
    
end


%%