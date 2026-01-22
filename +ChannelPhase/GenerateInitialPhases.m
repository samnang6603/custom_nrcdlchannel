function [Phi,cdl_struct] = GenerateInitialPhases(cdl_struct,ant_struct,info_struct)

randomStream = cdl_struct.RandomStreamObj;

c0 = 299792458; % Light speed

clusterTypes = info_struct.ClusterTypes;

M = RayModel.GetNumberOfRays;
L = length(clusterTypes);
LOS = find(strcmpi(clusterTypes,'LOS')); % index of LOS
splitNLOS = find(strcmpi(clusterTypes,'SubClusteredNLOS'));

% Phi Size: The 4 refers to 4 polarization combinations
% (theta-theta,theta-phi,phi-theta,phi-phi), theta: zenith, phi: azimuth
% See TR 38.901 section 7.5 for more details
Phi_size = [L M 4];

if (strcmpi(cdl_struct.InitialPhases,'Random'))
    cdl_struct.InitialPhases = '38.900';
end

if any(strcmpi(cdl_struct.InitialPhases,{'38.900','38.901','36.873'}))
    if isempty(randomStream)
        Phi = rand(Phi_size);
    else
        Phi = rand(randomStream,Phi_size);
    end
    % rand() is uniform from [0,1) but we need (-pi,pi)
    Phi = Phi*2*pi - pi;
elseif isscalar(cdl_struct.InitialPhases)
    Phi = repmat(double(cdl_struct.InitialPhases),Phi_size);
end

if ~isempty(LOS) % has LOS component
    Phi(LOS,2:M,:) = -Inf; % LOS only has 1 ray, realistically, so replaces other rays with -Inf
    Phi(LOS,1,2:4) = -Inf; % only has phase on the theta-theta combo first ray for LOS cluster
    if strcmpi(cdl_struct.InitialPhases,'38.901') || ...
            strcmpi(cdl_struct.InitialPhases,'38.900') || ...
            isnumeric(cdl_struct.InitialPhases) ...
        % See TR 38.901 Equation 7.5-29
        % Phase of exponential term with d_3D
        lambda_0 = c0/cdl_struct.CarrierFrequency; % wavelength of carrier tone

        % d_3D: effective projection of element separation in the ray's
        % direction of arrival/departure
        % Imagine the ray as a laser beam coming from a certain direction 
        % (described by zenith and azimuth angles). 
        % Then imagine your antenna element sitting somewhere in space.
        % Now:
        % Drop a perpendicular from that antenna element onto the laser 
        % beam's direction vector. Measure that distance along the 
        % ray direction. That's d_3D.
        % It tells you:
        % How far "along" the ray this element is located.
        % Which in turn tells you how much phase shift this particular ray 
        % will accumulate by the time it hits (or leaves) this element.
        % Thus d_3D is a 3D Euclidean distance between position of Tx and
        % Rx: d_3D = ||r_Rx - r_Tx|| norm-2
        d_3D = sqrt(sum((ant_struct.TransmitAntennaArray.Position -...
            ant_struct.ReceiveAntennaArray.Position).^2));
        
        % The phase
        Phi(LOS,1) = -2*pi*d_3D/lambda_0;

    end
end

end