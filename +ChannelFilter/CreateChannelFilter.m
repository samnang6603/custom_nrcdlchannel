function filt_struct = CreateChannelFilter(cdl_struct,sigin)

filt_struct = cdl_struct.ChannelFilter;
filt_struct = designFilterSpecs(filt_struct,cdl_struct,sigin);

end

function filt_struct = designFilterSpecs(filt_struct,cdl_struct,sigin)

% Sample Rate
sampleRate = filt_struct.SampleRate;

% Channel filter uses 16 FIR taps by default in comm.ChannelFilter.
% nrCDLChannel sets delay to 7 samples (i.e., floor of group delay 7.5).
% This balances minimal delay and causality while avoiding misalignment.
% The 16-tap FIR is fixed and optimized to smooth time-varying channel
% coefficients without introducing extra latency beyond what's necessary.
filtLen = 16;

% Interpolation factor 'L' is chosen based on desired fractional delay error.
% To bound worst-case error in delay interpolation to 0.01 samples,
% we use: L = ceil(1/(2*FracDelayError)), based on interpolation theory.
% This ensures tight temporal alignment between channel and waveform samples.
% Because in the worst-case, the interpolation error could land midway 
% between two interpolation points, that’s a half step in the fractional 
% grid.
L = ceil(1/(2*filt_struct.MaxFractionalDelayError));

% Path delays in samples
pathDelays = cdl_struct.PDP.PathDelays;
tpathDelays = pathDelays*sampleRate;
tpathDelays = tpathDelays(:);

% Fractional delays in samples
fracDelay = mod(tpathDelays,1);

% Interpolation phases ranging from 1 to L+1
% Define the fractional delay interpolation phases in reverse order.
% This matches MATLAB’s internal polyphase filter layout, where taps for
% fractional delays are stored starting from largest (near 1) to smallest.
% Ensures correct indexing into filter bank and tap alignment for interpolation.
interpRange = [0, (L-1:-1:1)/L, 1];
[~,phaseIdx] = min(abs(fracDelay - interpRange),[],2);

% If fracDelay is closest to 0 or 1, do not interpolate
pathOnSampleTime = any(phaseIdx == [1 L+1],2);

% Calculate integer delay part
intDelay = floor(tpathDelays) + (phaseIdx == L+1);

% Tap index range, excluding the path with fracDelay closest to 0 or 1
tapRange = intDelay + (~pathOnSampleTime.*[-filtLen/2+1,filtLen/2]);

% Tap index range to cover all the paths.
tapDelay = min([0;tapRange(:,1)]):max(tapRange(:,2));

% *** Info on Tap delay and channel filter delay ***
% - Tap delay refer to propagation environment (CDL). It affects actual
% multipath spacing. Each tap is one cluster/echo path with a certain delay
% Comes from PDP. Unit is usually nanosecond.
% - Channel filter delay refer to FIR filter (interpolation/smoothing). It
% affects only the filtering mechanics and output alignment. Unit is
% usually in samples.
% *************************************************************************

% Channel filter delay
channelFilterDelay = filt_struct.ChannelFilterDelay;

% Update integer delays, get the tstart from tapRange and add in the
% channel filter delay
pathIntDelays = tapRange(:,1) + channelFilterDelay;

% Configure stopband attenuation
stopbandAttenuation = filt_struct.StopbandAttenuation;
if ~all(pathOnSampleTime)
    stopbandAttenGrid = repmat(stopbandAttenuation, 1, L);
    fracDelayCoeff = designInterpFilter(L,filtLen,phaseIdx,stopbandAttenGrid);
else
    fracDelayCoeff = 1;
end

[filtCoeff_cell,filtState_cell] = getFilterCoeff(sigin, pathOnSampleTime,...
    fracDelayCoeff,pathIntDelays);

filt_struct.FilterCoefficient = filtCoeff_cell;
filt_struct.FilterState = filtState_cell;

end

function interpMatrix = designInterpFilter(L,filtLen,phaseIdx,stopbandAttenGrid)
nphaseIdx = length(phaseIdx);
interpMatrix = zeros(nphaseIdx,filtLen);
for m = 1:nphaseIdx
    if any(phaseIdx(m) == [1,L+1])
        interpMatrix(m,filtLen/2) = 1; % unused
    else
        polyphaseFIR = designMultirateFIR(L,1,filtLen/2,stopbandAttenGrid(phaseIdx(m)));
        h = reshape(polyphaseFIR,L,[]);
        interpMatrix(m,:) = h(phaseIdx(m),:);
    end
end
end

function [filtCoeff_cell,filtState_cell] = getFilterCoeff(sigin, pathOnSampleTime, fracDelayCoeff, intDelay)

% Create array of cell for easier access and operation
fracDelayCoeff_cell = mat2cell(fracDelayCoeff,...
    ones(1, size(fracDelayCoeff,1)), size(fracDelayCoeff,2));

% If path is on sample time, then just replace with 1
fracDelayCoeff_cell(pathOnSampleTime) = {1};

% Allocate filter coeff and state
filtCoeff_cell = cell(size(fracDelayCoeff_cell));
filtState_cell = filtCoeff_cell;

% Loop over every cell element and 
for m = 1:length(fracDelayCoeff_cell)
    %pad in zero the excess integer delay 
    % with the corresponding fracDelayCoeff when necessary
    tmp = [zeros(1, intDelay(m)), fracDelayCoeff_cell{m}];
    tmpLength = length(tmp);

    % Assign filter coeff and state
    filtCoeff_cell{m} = tmp;
    filtState_cell{m} = zeros(tmpLength-1,size(sigin,2),size(sigin,3),like=1i);
end


end