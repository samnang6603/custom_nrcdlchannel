function mdl_struct = CreatePowerDelayProfileStructure(cdl_struct)
% CREATEPOWERDELAYPROFILESTRUCTURE Extract PDP and scaled with some
% parameters, then extract per cluster parameter.
%
%

delaySpread = cdl_struct.DelaySpread;
delayProfile = cdl_struct.DelayProfile;

% Get CDL_A, CDL_B, CDL_C, CDL_D, CDL_E, LUTs
pdp_table = PowerDelayProfile.GetDelayProfileTable(delayProfile);

% For CDL-D, CDL-E, has LOS
desiredKFactor = NaN;
hasLOS = PowerDelayProfile.CheckForLoS(cdl_struct);
if hasLOS
    mdl_struct.HasLoSCluster = true;

    % If manually change K Factor value
    if cdl_struct.KFactorScaling 
        desiredKFactor = cdl_struct.KFactor;
    end
else
    mdl_struct.HasLoSCluster = false;
end

% Scale target K-Factor and delay spread 
pdp_table_out = PowerDelayProfile.ScaleDelaysKFactor(pdp_table,desiredKFactor,delaySpread);

% Assign to output struct 
% The column indices are in the CDL profile tables in TR 38.901 7.7.1
% The column indices are listed below
%     0     |         1        |   2   |  3  |  4  |  5  |  6
% Cluster # | Normalized Delay | Power | AoD | AoA | ZoD | ZoA 
% 
% For CDL-D/E, the Specular(LoS) component is in cluster row #1
mdl_struct.PathDelays = pdp_table_out(:,1).';
mdl_struct.AveragePathGains = pdp_table_out(:,2).';
mdl_struct.AnglesAoD = pdp_table_out(:,3).';
mdl_struct.AnglesAoA = pdp_table_out(:,4).';
mdl_struct.AnglesZoD = pdp_table_out(:,5).';
mdl_struct.AnglesZoA = pdp_table_out(:,6).';

% Overall table for completeness
mdl_struct.PDPTable = pdp_table_out;

% Assign other info to output struct
mdl_struct.NormalizedPathGains = cdl_struct.NormalizePathGains;
mdl_struct.DelaySpread = delaySpread;
mdl_struct.DelayProfile = delayProfile;

% Angles Scaling
mdl_struct.AngleScaling = cdl_struct.AngleScaling;
if (mdl_struct.AngleScaling)
    MeanAngles = cdl_struct.MeanAngles;
    mdl_struct.TargetMeanAoD = MeanAngles(1);
    mdl_struct.TargetMeanAoA = MeanAngles(2);
    mdl_struct.TargetMeanZoD = MeanAngles(3);
    mdl_struct.TargetMeanZoA = MeanAngles(4);
end

% Get per cluster parameter and assign to output struct
per_cluster = PowerDelayProfile.GetPerClusterParam(cdl_struct);
mdl_struct.XPR = per_cluster.XPR;
mdl_struct.AngleSpreads = [per_cluster.C_ASD,... 
                          per_cluster.C_ASA,...
                          per_cluster.C_ZSD,... 
                          per_cluster.C_ZSA];
mdl_struct.NumStrongestCluster = 0;


end