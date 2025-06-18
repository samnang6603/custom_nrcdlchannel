function F = ComputePolarization(polmode,powpattern_fcn_handle,powmode,theta_p,phi_p,zeta)
% POLARIZATIONMODE Select polmode based on Model-1 or 2 as a function
% of the implicit function powpattern_fcn_handle, which is A_prime, power
% radiation pattern and the angles
% Note: The '_p' denotes 'prime' status, '_pp' denotes 'double prime'
% Inputs:
%        - polmode: Polarization Model ('Model-1', 'Model-2')
%        - powpattern_fcn_handle: implicit A_prime, power radiation pattern
%        - powmode: Power pattern mode ('38.901'or 'BS', 'isotropic')
%        - theta_p: Zenith angle
%        - phi_p: Azimuth angle
%        - zeta: Polarization angles
% Output:
%        - F: Field pattern of element
switch lower(polmode)
    case 'model-1'

        % TR 38.901 7.3.2 Model-1
        % Pol field effect of Zenith component as a fcn of LCS angle
        F_theta_pp = sqrt(powpattern_fcn_handle(powmode,theta_p,phi_p));

        % Rotation matrix elements cos(Psi) and sin(Psi) for an angular
        % displacement of Psi due to the orientation of the LCS w.r.t. the GCS
        % See Equation 7.3-3 cos(phi) and sin(phi)
        denominator = sqrt(1 - ((cosd(zeta).*cosd(theta_p)) - (sind(zeta).*sind(phi_p).*sind(theta_p))).^2);
        cosPsi = ((cosd(zeta).*sind(theta_p)) + (sind(zeta).*sind(phi_p).*cosd(theta_p))) ./ denominator;
        sinPsi = (sind(zeta).*cosd(phi_p)) ./ denominator;

        % Assume vertical polarization in cases where the transformation
        % degenerates, this is, when theta_p is near zeta or 180-zeta and
        % phi_p is near -90 or 90, respectively. When zeta = 0/180, any
        % value of phi_p makes the transformation degenerate. The threshold
        % (10^-5) is the upper bound of the magnitude error outside a region of
        % radius ~10^-4 degrees around the singularity. Within that region, the
        % polarization angle error can be arbitrary.
        degen = isnan(cosPsi) | isnan(sinPsi) | ...
            ( abs(hypot(sinPsi,cosPsi) - 1) > 1e-5 );
        if any(degen(:)) % degenerates protection
            cosPsi(degen) = 1;
            sinPsi(degen) = 0;
        end

        % Equation 7.3-3 evaluated assuming F_phi_pp = 0
        %{
            [F_theta_p]  =  [cos(phi) -sin(phi)][F_theta_pp]
            [F_phi_p  ]     [sin(phi)  cos(phi)][F_phi_pp]
        %}
        F = [F_theta_pp.*cosPsi; F_theta_pp.*sinPsi];

    case 'model-2'
        sqrtA_p = sqrt(powpattern_fcn_handle(powmode,theta_p,phi_p));

        % TR 38.901 Equation 7.3-4
        F_theta_p = sqrtA_p*cosd(zeta);

        % TR 38.901 Equation 7.3-5
        F_phi_p = sqrtA_p*sind(zeta);

        % Now combine
        F = [F_theta_p; F_phi_p];

end
