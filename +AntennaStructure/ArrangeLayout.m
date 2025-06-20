function ant_layout = ArrangeLayout(arrayStruct,orientation)
% Inspired by MATLAB Antenna layout arrangement in nrCDLChannel()

% Initial array position
ant_layout.Position = [0; 0; 0];

% Initial array orientation
ant_layout.Orientation = orientation;

% Dimension sizes [M N P M_g N_g] of antenna array
%{
    M: Number of rows of antenna elements
    N: Number of columns of antenna elements
    P: Number of possible polarization config
    M_g: Number of rows of antenna panels
    N_g: Number of columns of antenna panels
%}
dims = arrayStruct.Size; % array dimension [M N P M_g N_g]
num_elem = length(dims);
if num_elem ~= 5
    error('Antenna array must be specified in [N N P M_g N_g] format')
end
M = dims(1);
N = dims(2);
P = dims(3);
M_g = dims(4);
N_g = dims(5);

% Antenna element and panel spacing 
elem_spacing = arrayStruct.ElementSpacing;
dV = elem_spacing(1); % vertical spacing between elements
dH = elem_spacing(2); % horizontal spacing between elements
dgV = elem_spacing(3); % vertical spacing between panels
dgH = elem_spacing(4); % horizontal spacing between panels

% element polarization zeta
zeta = arrayStruct.PolarizationAngles;

% Create a matrix that lay out the array dimensions
% [M x N x P x M_g x N_g]. The third (polarization) column vector has no
% spacing
s = [[dV; 0; 0], [0; dH; 0], [0; 0; 0], [dgV; 0; 0], [0; dgH; 0]];

% Create positions matrix for elements, pol, and panels
pos = cell(1,num_elem);
for m = 1:num_elem
    % repeat column vector here according to number of
    % antenna elements/polarizations/panels
    pos{m} = repmat(s(:,m),1,dims(m));

    % now map out local coordinate and center the array
    tmp = cumsum(pos{m},2);
    center = tmp - mean(tmp,2);
    pos{m} = center;
end

% Expand the positions for antenna elems, pol, and panels into higher 
% dimension array to account for all effects and interactions between:
% - element and element
% - element and polarization
% - cross-panel element interaction
% spm is called spatial map
spm = ExpandSpatialCorrTensor(pos,[3 M N P M_g N_g]);

% Combine all the elements in cell array spm to create one spatial map
% Stack 3D element tensors into one and reduce
spmstacked = cat(num_elem+2, spm{:}); % concat to expand to next new dim
spmsum = sum(spmstacked, num_elem+2);

% Define reorientation matrix for the array's local broadside direction
% By default, this array is defined to point toward the +x-axis in its own 
% LCS.3GPP, however, requires that the array broadside points toward the 
% +z-axis in the GCS. This matrix rotates the local +x direction to align 
% with global +z.
yz2xy_rot = [0 0 1; 0 1 0; -1 0 0];  % Rotation from local-x to global-z

% Build the orientation matrix from Euler angles
% Uses TR 38.901 composite rotation matrix R = Rz(bearing)*Ry(downtilt)*Rx(slant)
% Converts the LCS into GCS using active rotation.
vect_gcs = AntennaStructure.TransformVectorLCS2GCS(ant_layout.Orientation);  % 3x3 orientation matrix

% Combine reorientation and orientation matrix
% First rotate from local-x to global-z (yz2xy_rot)
% Then rotate according to the array's orientation in the scenario
vect_gcs_zbroadside = vect_gcs*yz2xy_rot;  % Final LCS → GCS mapping

% Apply the transformation to all stacked element positions
% Reshape the 6D tensor [3 x M x N x P x Mg x Ng] into a 2D matrix [3 x total_elements]
spmsum_mat = reshape(spmsum, 3, prod([M N P M_g N_g]));

% Rotate every position vector from LCS into its GCS equivalent
reorient = vect_gcs_zbroadside*spmsum_mat;

% Reshape back to the original 6D tensor format
reorient = reshape(reorient, [3 M N P M_g N_g]);

% Save transformed positions
ant_layout.ElementPositions = reorient;

% 'zeta' holds the polarization slant angles (e.g. [45 -45] for ± slant)
% We expand it across all element positions,
zeta_tensor = repmat(zeta',[1 M N M_g N_g]); % turn zeta in tensor

% Orientation affects polarization of antenna elements
% so that the polarization angle sits in the 3rd vector component of
% ormap (ormap is [3 x M x N x P x Mg x Ng])
ormap = zeros([3 M N P M_g N_g]);
ormap(3,:,:,:,:,:) = permute(zeta_tensor,[2 3 1 4 5]); % put zeta into respective 3rd dim for pol
% Why permute(2 3 1 4 5)? Because:
% - 'zeta_tensor' was [2 x M x N x Mg x Ng]
% - We want [M x N x 2 x Mg x Ng] to align it correctly under ormap(3,...)
ant_layout.ElementOrientations = ormap;

% Preallocate extra stuffs
ant_layout.SubarrayPositions = [];
ant_layout.InitialArrayOrientation = [0; 0; 0];
end

function y = ExpandSpatialCorrTensor(x,arraydims)
%EXPANDSPATIALCORRTENSOR Expands spatial position tensors for antenna 
% elements.
%
%   y = ExpandSpatialCorrTensor(x, arraydims) takes a cell array of spatial 
%   vectors x (typically 3x1 vectors for [y; x; z] positions) and expands 
%   them into higher-dimensional arrays that match the full antenna array 
%   configuration specified in arraydims = [M, N, P, M_g, N_g].
%
%   The function performs:
%       - Repmat expansion of each spatial axis into all other dimensions
%       - Permutes the resulting tensor so that its dimensions match the 
% original
%         array layout for correct spatial correlation calculations
%
%   Inputs:
%       x          - Cell array of spatial position vectors (e.g., x, y, z)
%       arraydims  - Antenna array dimensions in the form [3 M N P Mg Ng]
%                    Number 3 represents y,x,z local coordinate
%
%   Output:
%       y          - Cell array with each entry expanded to the shape 
%                    matching arraydims, suitable for spatial pattern or 
%                    correlation modeling
y = x;
currentdim_idx = 1;
for m = 1:length(x)
    repdims = arraydims;

    % Index of current spatial axis (e.g., y, x, or z)
    currentdim_idx = currentdim_idx + 1;
    active_dims = [1 currentdim_idx];
    repdims(active_dims) = []; 
    
    % Expand along remaining antenna dimensions
    y{m} = repmat(x{m}, [ones(1, length(active_dims)) repdims]);

    % Reorder dimensions to match original array layout
    permdims = 1:length(arraydims);
    permdims(active_dims) = [];
    permdims = [active_dims permdims]; %#ok<AGROW>
    [~, sortidx] = sort(permdims);
    y{m} = permute(y{m}, sortidx);
end

end

