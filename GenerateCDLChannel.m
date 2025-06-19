function mdl_out = GenerateCDLChannel(cdl_struct)
% Generate NR CDL Channel, inspired by MATLAB nrCDLChannel()

% Step 1 configure antenna structure
ant_struct = AntennaStructure.ConfigureStructure(cdl_struct);

% Step 2 configure delay profile
pdp_struct = PowerDelayProfile.CreatePowerDelayProfileStructure(cdl_struct);

% Step 3 split LoS an NLoS cluster and sub-clustering
[pdp_struct,infostruct] = PowerDelayProfile.GetLOSClusterInfo(pdp_struct);

%iniphase_struct;


mdl_out = 0;
end