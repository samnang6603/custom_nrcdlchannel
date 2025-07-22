function sigout = FilterCDLChannel(sigin,H,sampleTimes,cdl_struct)

% Get input and output sizes
insiz = size(sigin);
numTimeSamples = insiz(1);
numInputSignals = insiz(2);
numOutputSignals = size(H,4);

% Allocate output size
sigout = zeros(numTimeSamples,numOutputSignals,like=1i);

% Get pathGain indices
[idx,endpts] = pathGainIndices(cdl_struct.SampleDensity,sampleTimes,insiz,cdl_struct.SampleRate);

% Create filter structure with property
filt_struct = ChannelFilter.CreateChannelFilter(cdl_struct,sigin);

% Set filter mode for memory management
[filt_struct,endpts] = configureFilterMode(H,idx,endpts,filt_struct);

% Number of filter sections
nFilterSections = size(endpts,1);

if nFilterSections == 1 && (size(H,1) == numTimeSamples || size(H,1) == 1)
    % Mode (1)
    % If there is only one section and the path gain matches the input
    % sampling rate exactly/are scalar, then do channel filtering directly
    sigout = channelFilterSubroutine(filt_struct,cdl_struct,sigin,H);
else
    % By sections
    for m = 1:nFilterSections
        tstart = endpts(m,1);
        tend = endpts(m,2);
        if filt_struct.FilterPolicy.FilterMode == 1
            % If section but remains in Mode (1) (over 1GiB memory per
            % waveform section
            Hm = H(idx(tstart:tend),:,:,:);
        else
            % Mode (2)
            Hm = H(m,:,:,:);
        end
        sigout(tstart:tend,:) = channelFilterSubroutine(filt_struct,...
            cdl_struct,sigin(tstart:tend,1:numInputSignals),Hm);
    end
end

end

function [idx,endpts] = pathGainIndices(sampleDensity,sampleTimes,inputSize,sampleRate)

if sampleDensity == inf

    T = inputSize(1);
    idx = 1:T;
    idx = idx(:);
    endpts = [idx, idx];
else
    % -- Align input waveform samples to time-varying channel matrices --
    % The channel matrix H(t) is computed at coarse time resolution sampleTimes,
    % while the input waveform is sampled at fine resolution (t).
    % We apply Zero-Order Hold (ZOH), where each waveform sample t(i) uses
    % the most recent past channel sample sampleTimes(j), 
    % i.e., H(t(i)) = sampleTimes(j) for sampleTimes(j) <= t(i).
    % This block computes a mapping (idx) and endpoints for efficient 
    % block-wise filtering.
    [idx,endpts] = getZOHChannelMatrix(sampleTimes,inputSize(1),sampleRate);

end

end

function [idx,endpts] = getZOHChannelMatrix(HsampleTimes,inputSize,sampleRate)

% Get time points
t = (0:inputSize-1).'/sampleRate;

% Move the H sample time points to the center of two sample time points by
% using mean

LH = length(HsampleTimes);

if LH > 1
    HsampleTimes = HsampleTimes + mean(diff(HsampleTimes))/2;
end

% Inlude t(1) in the H sample time points
if ~isempty(t)
    if t(1) < HsampleTimes(1)
        HsampleTimes = [t(1); HsampleTimes];
    end
end

% Align each input waveform sample time t with a channel matrix
% HsampleTimes.
% Channel coefficients H are only computed at discrete times HsampleTimes 
% (low rate) Waveform samples occur at finer time resolution t (high rate)
% We perform Zero-Order Hold (ZOH) to associate each t(i) with its
% most recent channel coefficient sample in HsampleTimes.

if LH > 1
    % Use histogram binning to map each waveform sample time t(i) to the
    % bin defined by channel time points HsampleTimes
    % - idx(i) contains the channel matrix index to use for sample t(i)
    % - endpts(k,:) identifies contiguous sections of t that share the same
    % HsampleTimes
    [endpts,~,idx] = histcounts(t,HsampleTimes);

    % Remove empty bins, then compute start and end indices of each bin group
    endpts(endpts==0) = [];
    endpts = [cumsum([1; endpts(1:end-1).']), cumsum(endpts.')];
else
    % If only a single channel sample time exists, apply it across all input samples
    idx = ones(inputSize,1);
    if ~inputSize
        endpts = zeros(0,2);
    else
        endpts = [1, inputSize];
    end
end

end

function [filt_struct,endpts] = configureFilterMode(H,idx,endpts,filt_struct)
% Ncs is number of Channel Snapshot or waveform section
% Np is the number of path
% NTx is the number of transmitter or input signals
% NRx is the number of receiver or output signals
[Ncs,Np,NTx,NRx] = size(H); 
bytesPerElement = 16; % only double or 64-bit floating pt precision

% Memory size that will be allocated for H after ZOH
memreq = prod([length(idx), Np, NTx, NRx])*bytesPerElement;

% memory per waveform section used by interpolated H
memsec = memreq/Ncs;

% Choose between two channel filtering modes:
% Mode (1): ZOH-interpolate full g(t) -> 1 filter call -> faster, high RAM
% When the memsec size exceeds 1 GiB, split into sections but still 
% retain Mode (1) processing
% Mode (2): Filter section-by-section without full ZOH -> lower RAM, more calls
%
% We switch to Mode (2) when full interpolation exceeds ~1 GiB
% to avoid memory overflows in large-scale MIMO or long waveform sims.

if Ncs == 1
    filterMode = 2;
else
    % According to MATLAB benchmark, memsec = 3 MiB seems to be the 
    % threshold where filter mode (1) and (2) are equally viable. Anything 
    % larger than this, (2) is preferred and anything lower, (1) is 
    % preferable.
    filterMode = (memsec > 3*2^20) + 1;
end

% Mode (1) but split when blocksize exceed 1 GiB
if filterMode == 1 && ~isempty(idx)
    num1GiB = ceil(memreq/2^30);
    etmp = floor(linspace(0,numel(idx),num1GiB+1));
    etmp = etmp(:);
    endpts = [etmp(1:end-1)+1, etmp(2:end)];
end

filt_struct.FilterPolicy.Ncs = Ncs;
filt_struct.FilterPolicy.MemoryRequired = memreq;
filt_struct.FilterPolicy.MemoryPerWaveformSection = memsec;
filt_struct.FilterPolicy.FilterMode = filterMode;

end

function y = channelFilterSubroutine(filt_struct,cdl_struct,x,h)
Ns = size(x,1);
NTx = size(x,2);
NRx = size(h,4);
Nb = size(x,3);
Npath = length(cdl_struct.PDP.PathDelays);

% Allocate
y = zeros(Ns,NRx,Nb,like=x);

% Get filter coeff and state
filtCoeff_cell = filt_struct.FilterCoefficient;
filtState_cell = filt_struct.FilterState;


if size(h,1) ~= 1 % Filter Mode (1) split, may remove later if never used
    for q = 1:Npath
        % This is a polyphase FIR from designInterpFIR from 
        % CreateChannelFilter
        b = filtCoeff_cell{q}; 
        % Apply polyphase interpolation filtering, yy is the filtered 
        % output or interpolated samples
        [yy, filtState_cell{q}] = filter(b,1,x,filtState_cell{q});
        
        % Apply path gain to yy
        for u = 1:NRx
            tmp = y(:,u,:);
            for s = 1:NTx
                hqus = reshape(h(:,q,u,s,:),Ns,1,Nb);
                yy2 = yy(:,s,:).*hqus;
                tmp = tmp + yy2;
            end
            y(:,u,:) = tmp;
        end    
    end
else
    % Channel gain
    hGain = permute(h,[3 4 2 1]); % Bring NTx to the front

    % Interpolate
    for q = 1:Npath
        % This is a polyphase FIR from designInterpFIR from
        % CreateChannelFilter
        b = filtCoeff_cell{q};
        % Apply polyphase interpolation filtering, yy is the filtered
        % output or interpolated samples
        [yy, filtState_cell{q}] = filter(b,1,x,filtState_cell{q});
        y = y + yy*hGain(:,:,q);
    end
end

end