function [H, sampleTimes] = GenerateTimeVaryingCDLChannel(cdl_struct,ant_struct,static_struct,scattererStates,scattererSpeeds)
%GENERATETIMEVARYINGCDLCHANNEL Add in Doppler term to static CDL channel

% Get time samples
sampsize = cdl_struct.NumTimesSamples;
P = ant_struct.NumInputSignals;
insize = [sampsize, P];
sampleTimes = getSampleTimes(cdl_struct,insize);

T = length(sampleTimes);
R = ant_struct.NumOutputSignals;
D = [T 1 P R];

staticPathGains = static_struct.StaticPathGains;
N = size(staticPathGains,2);
D(2) = N;

% Allocate channel matrix
H = zeros(D,like=1i);

% Get Doppler term according to TR 38.901 Equation 7.5-28/29 in step 11
txSphericalUnitVectors = static_struct.txSphericalUnitVectors;
rxSphericalUnitVectors = static_struct.rxSphericalUnitVectors;
doppler = computeDopplerTerm(cdl_struct,rxSphericalUnitVectors,...
    txSphericalUnitVectors,scattererStates,scattererSpeeds,sampleTimes);

% Apply doppler effect of each ray in the static path gains by multiplying
% the time-varying ray complex gains for each cluster

for u = 1:D(4) % Rx
    A = doppler;
    B = staticPathGains(:,:,:,:,u);
    % A and B shows different dims
    % A is size [#SampleTime x #Cluster x #Ray]
    % B is size [1 x #Cluster x #Ray x #InputSignals]
    % C = A.*B is size [#SampleTime x #Cluster x #Ray x #InputSignals]
    % For Example, let
    % A be size [2 x 14 x 20], B be [1 x 14 x 20 x 8]
    % MATLAB implicitly broadcasts singleton dims in A and B:
    %   - Expands dim 1 of B from 1 -> 2
    %   - Adds new 4th dim to A from 1 -> 8
    % Result: C is [2 x 14 x 20 x 8]
    C = A.*B;
    sumTerm = sum(C,3); % 3rd dim is the cluster dim
    % permute to per ray dim
    % H is [#SampleTime x #Cluster x #InputSignals x #Ray]
    H(:,:,:,u) = permute(sumTerm,[1 2 4 3]); 
end



end

function [timeSamp,F_cg] = getSampleTimes(cdl_struct,insize)
%GETSAMPLETIMES Return time samples based on cdl_struct param and
% cofficient generation sampling rate F_cg

c0 = 299792458;
fc = cdl_struct.CarrierFrequency;
initTime = cdl_struct.InitialTime;
maxDoppler_v = cdl_struct.MaximumDopplerShift; % Hz
if isscalar(maxDoppler_v)
    maxDoppler_v(2) = 0;
    maxScatter_v = 0; % depends on MaximumDopplerShift
else
    maxScatter_v = cdl_struct.MaximumScattererSpeed;
end

inputSampRate = cdl_struct.SampleRate;
sampDensity  = cdl_struct.SampleDensity;
numWaveformSamples = insize(1);

if isfield(cdl_struct,'SampleTimes')

    % Use sample times as time points
    timeSamp = cdl_struct.SampleTimes(:);
    F_cg = [];

else

    % Set F_cg according to Doppler frequency and SampleDensity param
    lambda_0 = c0/fc;

    % Sum total Doppler shift from both tx and rx
    txrxDoppler = sum(maxDoppler_v); % Hz

    % Doppler from moving scatterer, causing both the forward and reflected
    % path, like a double bounce
    scattererDoppler = 2*maxScatter_v/lambda_0; % Hz

    % The number 2 indicates F_cg must sample at least 2 times the max
    % Doppler freq as per Nyquist theorem but can be oversampled by
    % sampDensity factor
    oversample_factor = 2*sampDensity;

    % Calculate coefficient generation sampling rate
    F_cg = oversample_factor*(txrxDoppler + scattererDoppler);

    % Calculate total waveform duration in seconds
    timeWaveform_sec = numWaveformSamples/inputSampRate;

    % Compute the number of channel coefficient samples needed to cover
    % the waveform duration, sampled at F_cg (channel generation rate)
    numTimePoints = timeWaveform_sec*F_cg;

    % Pad by at least half a sample period to ensure ZOH has one final
    % sample to hold — prevents end-of-waveform interpolation artifacts
    numTimePoints = ceil(numTimePoints + 0.5);

    % Oversample channel taps to resolve Doppler effects (2*fD*SampleDensity),
    % but clip to SampleRate since that's the max resolution the waveform
    % allows. Physics first, practical constraints second.
    F_cg = min(F_cg,inputSampRate);
    numTimePoints = min(numTimePoints,numWaveformSamples);

    % calculate time points at F_cg
    if numTimePoints == 1 % if zero Doppler (numTimePoints=1, timeSamp=0, F_cg=0)
        timeSamp = 0;
    else
        timeSamp = (0:(numTimePoints-1)).'/F_cg;
    end

    timeSamp = timeSamp + initTime;

end



end

function doppler = computeDopplerTerm(cdl_struct,rhat_rx,rhat_tx,scattererStates,scattererSpeeds,t)
% Compute Doppler term according to TR 38.901 Equation 7.5-28/29 in step 11
%{
                        T
           /        rhat_(rx)_(LOS/NLOS) * vbar    \
        exp|j2*pi--------------------------------*t|
           \                 lambda_0              /

%}
c0 = 299792458;
lambda_0 = c0/cdl_struct.CarrierFrequency;

% Reshape rhat_rx, rhat_tx, alph and D to combine cluster and ray dimensions into a single
% row
sphVectSiz = size(rhat_rx);
rhat_tx = reshape(rhat_tx,[3, numel(rhat_tx)/3]);
rhat_rx = reshape(rhat_rx,[3, numel(rhat_rx)/3]);
scattererStates = scattererStates(:); % Moving scatterer states
scattererSpeeds = scattererSpeeds(:);

% Rx and Tx speed
maxDopplerShift = cdl_struct.MaximumDopplerShift;
if isscalar(maxDopplerShift)
    maxDopplerShift(2) = 0;
end

v_rx = maxDopplerShift(1)*lambda_0; % rx speed
v_tx = maxDopplerShift(2)*lambda_0; % tx speed

% User-Terminal Direction of travel
[rUT,cUT] = size(cdl_struct.UTDirectionOfTravel);
if rUT == 2
    UTDirectionOfTravel = repmat(cdl_struct.UTDirectionOfTravel,1,2);
elseif cUT == 2
    UTDirectionOfTravel = repmat(cdl_struct.UTDirectionOfTravel',1,2);
end
theta_v_rx = UTDirectionOfTravel(2,1); % Rx zenith angle of travel
phi_v_rx = UTDirectionOfTravel(1,1);   % Rx azimuth angle of travel
theta_v_tx = UTDirectionOfTravel(2,2); % Tx zenith angle of travel
phi_v_tx = UTDirectionOfTravel(1,2);   % Tx azimuth angle of travel

v_rxSphericalUnitVector = AntennaStructure.GetSphericalUnitVector(phi_v_rx,theta_v_rx);
v_txSphericalUnitVector = AntennaStructure.GetSphericalUnitVector(phi_v_tx,theta_v_tx);

vbar_rx = v_rx*v_rxSphericalUnitVector;
vbar_tx = v_tx*v_txSphericalUnitVector;

% Compute doppler
vnm_rx = rhat_rx.'*vbar_rx;
vnm_tx = rhat_tx.'*vbar_tx;
scatterer = 2*scattererSpeeds.*scattererStates; % 2 because Tx->Scatterer->Rx
tmp = (vnm_rx + vnm_tx + scatterer)/lambda_0;

timeVaryingGrid = kron(t,(1i*2*pi*tmp).');
doppler = exp(timeVaryingGrid);
doppler = reshape(doppler,[length(t), sphVectSiz(2:end)]);

end
