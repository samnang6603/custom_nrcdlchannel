function coupling = ComputeRayCoupling(cdl_struct,clusterTypes)
%COMPUTERAYCOUPLING Calculate ray coupling according to TR 38.901 section
%7.5 step 8

randomStream = cdl_struct.RandomStreamObj;

M = RayModel.GetNumberOfRays;
L = length(clusterTypes);
LOS = find(strcmpi(clusterTypes,'LOS')); % index of LOS
splitNLOS = find(strcmpi(clusterTypes,'SubClusteredNLOS'));

% Size
% 3 because
% 1: AoD -> AoA coupling
% 2: AoD -> ZoA coupling
% 3. AoD -> ZoD coupling
coupling_size = [L M 3];

% Randomly generate coupling
if strcmpi(cdl_struct.RayCoupling,'Random')
    if isempty(randomStream)
        tmp = rand(coupling_size);
    else
        tmp = rand(randomStream,coupling_size);
    end
    [~,couplingTmp] = sort(tmp,2); % Sort weakest to strongest coupling within ray
    
    % Rearranging ray coupling to preserve directional consistency:
    % -------------------------------------------------------------
    % While ray angles (AoD, AoA, ZoD, ZoA) are randomly generated, 
    % we want to avoid complete chaos. 
    % Specifically, we preserve the sorted ZoD ↔ ZoA relationship across 
    % rays for each cluster.
    %
    % This maintains a logical pairing between departure and arrival 
    % directions (e.g., rays leaving high tend to arrive high), which 
    % reflects realistic scattering behavior.
    %
    % Without this, the model would assign arrival angles that have no 
    % directional correlation with their departure, breaking spatial 
    % symmetry, phase continuity, and any realistic MIMO behavior.
    %
    % NOTE: Even if RayCoupling = 'Random', this logic still maintains 
    % partial directional consistency by sorting ZoD -> ZoA mappings.
    coupling = zeros(coupling_size);

    % copy only AoA->AoD & AoD->ZoD coupling
    coupling(:,:,[1 3]) = couplingTmp(:,:,[1 3]);
    for m = 1:L % for all cluster
        % Here sort the AoA->ZoD in ascending order and get only the index
        [~,sortIdx] = sort(couplingTmp(m,:,3)); 

        % Map the sort index to rearrange AoA->ZoA according to AoA->ZoD
        coupling(m,:,2) = couplingTmp(m,sortIdx,2); 
    end

    % To be implemented later, for subcluster
end

% 
clusterTypeToAdjust = {'LOS';'SubclusteredNLOS';'NLOS'};
for m = 1:length(clusterTypeToAdjust)
    idx = strcmpi(clusterType,clusterTypeToAdjust{m});

    % Total number of cluster of the same type;
    numClustersOfType = sum(idx); 

    % Each cluster has local indices [1:M] for rays.
    % We now remap those to global indices across all clusters of this type
    % Think of each column as a ray index inside its own cluster.
    % To globalize it, we stride by cluster. That means:
    % For each ray index i in cluster k:
    % global_index = (i - 1) * stride + k
    stride = numClustersOfType;

    % Convert local coupling index into global indices
    stride_part = (coupling(idx,:,:)-1)*stride; % col-wise 

    % Offset to get correct cluster-wise position
    cluster_offset = repmat((1:numClustersOfType).',1,M,3);

    % Final remapped coupling matrix (global indices)
    coupling(idx,:,:) = stride_part + cluster_offset;

end


