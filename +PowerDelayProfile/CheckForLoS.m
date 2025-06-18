function hasLoS = CheckForLoS(cdl)
% CHECKFORLOS Check for LoS or Specular component
switch upper(cdl.DelayProfile)
    case {'CDL-A','CDL-B','CDL-C','TDL-A','TDL-B','TDL-C','NTN-TDL-A','NTN-TDL-B','NTN-TDLA100'}
        hasLoS = false;
    case {'CDL-D','CDL-E','TDL-D','TDL-E','TDLD30','TDLD10','NTN-TDL-C','NTN-TDL-D','NTN-TDLC5'}
        hasLoS = true;
    otherwise
        hasLoS = false;
end

end