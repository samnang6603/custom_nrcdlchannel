function mdl_out = GenerateCDLChannel(cdl)
% Generate NR CDL Channel, inspired by MATLAB nrCDLChannel()

% Step 1 configure antenna structure
ant_struct = AntennaStructure.ConfigureStructure(cdl);

% Step 2 configure delay profile
pdp_struct = PowerDelayProfile.CreatePowerDelayProfileStructure(cdl);

% Step 3 split LoS an NLoS cluster and sub-clustering
%pdp_info = 

%iniphase_struct;


mdl_out = 0;
end