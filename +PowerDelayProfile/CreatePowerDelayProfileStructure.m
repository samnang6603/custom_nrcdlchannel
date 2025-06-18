function mdlstruct = CreatePowerDelayProfileStructure(cdl)
% CREATEPOWERDELAYPROFILESTRUCTURE Extract PDP and scaled with some
% parameters, then extract per cluster parameter.
%
%

delaySpread = cdl.DelaySpread;
delayProfile = cdl.DelayProfile;

% Get CDL_A, CDL_B, CDL_C, CDL_D, CDL_E, LUTs
pdp_table = PowerDelayProfile.GetDelayProfileTable(delayProfile);

% For CDL-D, CDL-E has LoS
targetKFactor = NaN;
hasLoS = PowerDelayProfile.CheckForLoS(cdl);
if hasLoS
    targetKFactor = cdl.KFactor;
    mdlstruct.HasLoSCluster = true; 
else
    mdlstruct.HasLoSCluster = false;
end

% Scale target K-Factor and delay spread 
pdp = PowerDelayProfile.ScaleDelaysKFactor(pdp_table,targetKFactor,delaySpread);

% Assign to output struct 
% The column indices are in the CDL profile tables in TR 38.901 7.7.1
% The column indices are listed below
%     0     |         1        |   2   |  3  |  4  |  5  |  6
% Cluster # | Normalized Delay | Power | AoD | AoA | ZoD | ZoA 
% 
% For CDL-D/E, the Specular(LoS) component is in cluster row #1
mdlstruct.PathDelays = pdp(:,1).';
mdlstruct.AveragePathGains = pdp(:,2).';
mdlstruct.AnglesAoD = pdp(:,3).';
mdlstruct.AnglesAoA = pdp(:,4).';
mdlstruct.AnglesZoD = pdp(:,5).';
mdlstruct.AnglesZoA = pdp(:,6).';

% Assign other info to output struct
mdlstruct.NormalizedPathGains = cdl.NormalizePathGains;
mdlstruct.DelaySpread = delaySpread;
mdlstruct.DelayProfile = delayProfile;

% Angles Scaling
mdlstruct.AngleScaling = cdl.AngleScaling;
if (mdlstruct.AngleScaling)
    MeanAngles = cdl.MeanAngles;
    mdlstruct.TargetMeanAoD = MeanAngles(1);
    mdlstruct.TargetMeanAoA = MeanAngles(2);
    mdlstruct.TargetMeanZoD = MeanAngles(3);
    mdlstruct.TargetMeanZoA = MeanAngles(4);
end

% Get per cluster parameter and assign to output struct
per_cluster = PowerDelayProfile.GetPerClusterParam(cdl);
mdlstruct.XPR = per_cluster.XPR;
mdlstruct.AngleSpreads = [per_cluster.C_ASD,... 
                          per_cluster.C_ASA,...
                          per_cluster.C_ZSD,... 
                          per_cluster.C_ZSA];
mdlstruct.NumStrongestCluster = 0;


end