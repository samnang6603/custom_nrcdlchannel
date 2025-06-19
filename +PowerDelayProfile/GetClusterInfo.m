function [pdp_struct,infostruct] = GetClusterInfo(pdp_struct)
% Get LOS, NLOS and NLOS Subcluster informations


[LOS_pdp,K_dB] = getLOSClusterPathGains(pdp_struct);
NLOS_pdp = getNLOSClusterPathGains(pdp_struct);
NLOS_split_pdp = [];
clusterTypes = getClusterTypes(LOS_pdp,NLOS_pdp,NLOS_split_pdp);

pdpnew = [LOS_pdp; NLOS_split_pdp; NLOS_pdp];

% Update PDP
pdp_struct.PathDelays = pdpnew(:,1).';
pdp_struct.AveragePathGains = pdpnew(:,2).';
pdp_struct.AnglesAoD = pdpnew(:,3).';
pdp_struct.AnglesAoA = pdpnew(:,4).';
pdp_struct.AnglesZoD = pdpnew(:,5).';
pdp_struct.AnglesZoA = pdpnew(:,6).';
pdp_struct.PDPTable = pdpnew;

% Update info structure
infostruct.LOS_pdp = LOS_pdp;
infostruct.K_dB = K_dB;
infostruct.NLOS_pdp = NLOS_pdp;
infostruct.ClusterTypes = clusterTypes;

% Update info structure on angle spreads and XPR
infostruct.XPR = pdp_struct.XPR;
infostruct.AngleSpreads = pdp_struct.AngleSpreads;


end


function [LOS_pdp,K_dB] = getLOSClusterPathGains(pdp_struct)
% Extract LOS cluster and K Factor

% Get pdp table
pdptmp = pdp_struct.PDPTable;
powIdx = 2;

% Check for available LoS
LOS = pdp_struct.HasLoSCluster;
if LOS
    K_dB = pdptmp(1,powIdx) - pdptmp(2,powIdx);
else
    K_dB = -Inf;
end
% Do this indexing so that when no LOS, it produces empty vector
LOS_pdp = pdptmp((2-LOS):1,:); 
end

function NLOS_pdp = getNLOSClusterPathGains(pdp_struct)
% Extract NLOS clusters

% Get pdp table
pdptmp = pdp_struct.PDPTable;
powIdx = 2;

% To be completed later, need number of strongest cluster

NLOS_pdp = pdptmp(2:end,:);
end

function NLOS_split_pow = getSplitNLOSClusterPathGains(pdp_struct)
% To be completed later, need number of strongest cluster

end

function clusterTypes = getClusterTypes(LOS_pdp,NLOS_pdp,NLOS_split_pdp)
nLOS = size(LOS_pdp,1);
nSplitLOS = size(NLOS_split_pdp,1);
nNLOS = size(NLOS_pdp,1);

% Allocate cell for cluster types
clusterTypes = cell(nLOS + nSplitLOS + nNLOS,1);

for m = 1:length(clusterTypes)
    if m <= nLOS
        clusterTypes{m} = 'LOS';
    elseif m <= (nLOS + nSplitLOS)
        clusterTypes{m} = 'SubclusteredNLOS';
    else
        clusterTypes{m} = 'NLOS';
    end
end

end



