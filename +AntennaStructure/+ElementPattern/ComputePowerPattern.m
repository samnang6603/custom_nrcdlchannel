function A = ComputePowerPattern(powmode,theta_p,phi_p)
% COMPUTEPOWERPATTERN Compute A_prime, antenna gain pattern
% Antenna element radiation pattern is described in TR 38.901 Section 7.3.
% Inputs:
%        - powmode: Power pattern mode ('38.901'or 'BS', 'isotropic')
%        - theta_p: Zenith angle
%        - phi_p: Azimuth angle
% Outputs:
%        - A: power radiation pattern in linear scale

switch lower(powmode)
    case {'38.901','bs'}
        SLA_V = 30.0;     % Side-Lobe Attenuation (dB)
        theta_3dB = 65.0; % Vertical Half-Power Beam-Width (degrees)
        A_m = 30.0;       % Front-to-back attenuation ratio (dB)
        phi_3dB = 65.0;   % Horizontal Half-Power Beam-Width (degrees)
        G_max = 8.0;      % Maximum directional gain of an element (dBi)

        % antenna element vertical radiation pattern (dB)
        A_EV = -min(12*((theta_p-90)/theta_3dB).^2,SLA_V);

        % antenna element horizontal radiation pattern (dB)
        A_EH = -min(12*(phi_p/phi_3dB).^2,A_m);

        % combining method for 3D antenna element pattern (dB)
        A = -min(-(A_EV + A_EH),A_m);

        % incorporate maximum gain and convert to linear power
        A = 10.^((A + G_max)/10);

    case 'isotropic'
        % Equally radiate power in all directions, pattern is circle
        A = ones(size(theta_p));
end

end