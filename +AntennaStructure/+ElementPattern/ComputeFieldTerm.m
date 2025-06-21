function fieldTerm = ComputeFieldTerm(array_struct,theta,phi,antIdx)
%COMPUTEFIELDTERM Compute field term in GCS
% Using Equation 7.1-11. The process of this function involves finding psi
% and F_prime from polarization and field pattern
%

% Collapse and serialize theta and phi into single dim
theta = theta(:).';
phi = phi(:).';

% NOTE in TR 38.901:  
% - theta and phi are Azimuth and Zenith, respectively, in GCS
% - theta_prime and phi_prime as Azimuth and Zenith, respectively, in LCS
% This distinction is critical for transforming field vectors and applying 
% element pattern gains.

% Orientation angles of the array
% Bearing and Downtilt angles are the combination of array and element
% orientation
% Bearing angle (x-axis)
alphaInit = array_struct.InitialArrayOrientation(1); % initial bearing
alpha = array_struct.Orientation(1) + array_struct.ElementOrientations(1,antIdx);

% Downtilt angle (z-axis)
betaInit = array_struct.InitialArrayOrientation(2); % initial downtilt
beta = array_struct.Orientation(2) + array_struct.ElementOrientations(2,antIdx);

% Slant angle (y-axis)
gammaInit = array_struct.InitialArrayOrientation(3);
gamma = array_struct.Orientation(3);

% Rotate from LCS to GCS
orientInit = [alphaInit; betaInit; gammaInit];    % Initial/reference orientation (bearing, downtilt, slant)
orient     = [alpha; beta; gamma];               % Target/combined orientation (e.g., element + array)

% Compute rotation matrix from LCS to GCS for the *initial* orientation.
% NOTE: LCS2GCSRotationMatrix() returns R', so we take the transpose
% to obtain the proper LCS-to-GCS rotation matrix.
RInit = AntennaStructure.LCS2GCSRotationMatrix(orientInit)';

% Compute LCS-to-GCS rotation matrix for the *target* orientation.
ROriented = AntennaStructure.LCS2GCSRotationMatrix(orient);

% Transform the target orientation into the local frame of the initial 
% orientation.
% That is: R = RInit' * ROriented * RInit
% This performs a change of basis:
% - Maps ROriented (in GCS) into the LCS of the reference orientation.
% - The result, R, expresses the target orientation *as seen from* the 
%   reference LCS.
%
% Extra: To transform a matrix A into another basis B, we use
% A_in_B = B^T * A * B
% 
% WHY go back LCS frame?
% BECAUSE: We care about how the rotated antenna (element + array) looks 
% from the antenna’s starting orientation — because field vectors, 
% element patterns, and mutual coupling are defined in the antenna's own 
% local frame, not the world frame.
R = RInit' * ROriented * RInit;
Rp = R'; % For convenience for later computations, inverse = transpose

% Calculate rho, the unit radius in the spherical coordinate
% (rho,theta,phi) <=> (radius,azimuth,zenith)
% according to TR 38.901 Equation 7.1-6
sintheta = sind(theta); % for optimization purpose, reusable
costheta = cosd(theta);
rhohat = [sintheta.*cosd(phi)
          sintheta.*sind(phi)
               costheta     ];

% Compute theta_prime and phi_prime according to TR 38.901 Equation 7.1-7
% and 7.1-8
theta_p = real(acosd([0,0,1]*Rp*rhohat)); % azimuth in LCS
phi_p = atan2d(Rp(2,:)*rhohat,Rp(1,:)*rhohat); % phi in LCS
phi_p(theta_p==0) = 0; % Set ambiguous phi to 0 

% Next need to compute psi, the angular displacement between two pairs of
% unit vectors, according to Equation 7.1-12
% First, we need to find the unit vector of theta, phi, and theta_p. 
% These are the Cartesian representation of the spherical unit vectors

% Reusable constants for optimization
sinphi = sind(phi);
cosphi = cosd(phi);
costheta_p = cosd(theta_p);

% The unit vectors thetahat, phihat, and theta_phat
thetahat = [costheta.*cosphi
            costheta.*sinphi
               -sintheta    ];
phihat   = [-sinphi
             cosphi
             zeros(size(phi))];

theta_phat = [costheta_p.*cosd(phi_p)
              costheta_p.*sind(phi_p)
                  -sind(theta_p)     ];

% Using Equation 7.1-12, simplified by hand to get the least computational
% expression
psi = atan2d(sum(phihat.*(R*theta_phat)),sum(thetahat.*(R*theta_phat)));

% TR 38.901 Section 7.3.2
% Polarization slant angle
zeta = array_struct.ElementOrientations(3,antIdx);

% Get the Polarization and Element field pattern equation
polmode = array_struct.PolarizationModel;
powmode = array_struct.Element; 
F = array_struct.ElementPattern(polmode,powmode,theta_p,phi_p,zeta);

% Finally, we have everything we need to compute fieldTerm in GCS,
% using equation 7.1-11
cospsi = cosd(psi);
sinpsi = sind(psi);
%fieldTerm1 = [cospsi, -sinpsi; sinpsi cospsi]*F;

% Alternatively, simplied form, hand derived for optimization:

fieldTerm = [F(1,:).*cospsi - F(2,:).*sinpsi
             F(1,:).*sinpsi + F(2,:).*cospsi];
end