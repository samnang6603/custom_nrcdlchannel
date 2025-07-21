function pdp_struct = CreatePowerDelayProfileStructure(cdl_struct)
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
hasLOS = PowerDelayProfile.CheckForLOS(cdl_struct);
if hasLOS
    pdp_struct.HasLOSCluster = true;

    % If manually change K Factor value
    if cdl_struct.KFactorScaling 
        desiredKFactor = cdl_struct.KFactor;
    end
else
    pdp_struct.HasLOSCluster = false;
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
pdp_struct.PathDelays = pdp_table_out(:,1).';
pdp_struct.AveragePathGains = pdp_table_out(:,2).';
pdp_struct.AnglesAoD = pdp_table_out(:,3).';
pdp_struct.AnglesAoA = pdp_table_out(:,4).';
pdp_struct.AnglesZoD = pdp_table_out(:,5).';
pdp_struct.AnglesZoA = pdp_table_out(:,6).';

% Overall table for completeness
pdp_struct.PDPTable = pdp_table_out;

% Assign other info to output struct
pdp_struct.NormalizedPathGains = cdl_struct.NormalizePathGains;
pdp_struct.DelaySpread = delaySpread;
pdp_struct.DelayProfile = delayProfile;

% Angles Scaling
pdp_struct.AngleScaling = cdl_struct.AngleScaling;
if (pdp_struct.AngleScaling)
    MeanAngles = cdl_struct.MeanAngles;
    pdp_struct.TargetMeanAoD = MeanAngles(1);
    pdp_struct.TargetMeanAoA = MeanAngles(2);
    pdp_struct.TargetMeanZoD = MeanAngles(3);
    pdp_struct.TargetMeanZoA = MeanAngles(4);
end

% Get per cluster parameter and assign to output struct
per_cluster = PowerDelayProfile.GetPerClusterParam(cdl_struct);
pdp_struct.XPR = per_cluster.XPR;
pdp_struct.AngleSpreads = [per_cluster.C_ASD,... 
                          per_cluster.C_ASA,...
                          per_cluster.C_ZSD,... 
                          per_cluster.C_ZSA];
pdp_struct.NumStrongestCluster = 0;


end