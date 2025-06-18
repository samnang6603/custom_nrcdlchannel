function per_cluster = GetPerClusterParam(cdl)
% GETPERCLUSTERPARAM
% TR 38.901 Section 7.7.1 under each CDL table
switch upper(cdl.DelayProfile)
    case 'CDL-A'
        [C_ASD,C_ASA,C_ZSD,C_ZSA,XPR] = deal(5, 11, 3, 3, 10);
    case 'CDL-B'
        [C_ASD,C_ASA,C_ZSD,C_ZSA,XPR] = deal(10, 22, 3, 7, 10);
    case 'CDL-C'
        [C_ASD,C_ASA,C_ZSD,C_ZSA,XPR] = deal(2, 15, 3, 7, 7);
    case 'CDL-D'
        [C_ASD,C_ASA,C_ZSD,C_ZSA,XPR] = deal(5, 8, 3, 3, 11);
    case 'CDL-E'
        [C_ASD,C_ASA,C_ZSD,C_ZSA,XPR] = deal(5, 11, 3, 7, 8);
end
per_cluster.C_ASD = C_ASD;
per_cluster.C_ASA = C_ASA;
per_cluster.C_ZSD = C_ZSD;
per_cluster.C_ZSA = C_ZSA;
per_cluster.XPR = XPR;

end

