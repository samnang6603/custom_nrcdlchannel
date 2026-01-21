function [cdl_struct,ant_struct] = ConfigureStructure(cdl_struct)

TransmitAntennaArray = cdl_struct.TransmitAntennaArray;
ReceiveAntennaArray = cdl_struct.ReceiveAntennaArray;
TransmitArrayOrientation = cdl_struct.TransmitArrayOrientation;
ReceiveArrayOrientation = cdl_struct.ReceiveArrayOrientation;
CarrierFrequency = cdl_struct.CarrierFrequency;

mdlArrays = AntennaStructure.ArrangeAntennaStructure(...
    TransmitAntennaArray,...
    ReceiveAntennaArray,...
    TransmitArrayOrientation,...
    ReceiveArrayOrientation,...
    CarrierFrequency);

ant_struct = mdlArrays;
ant_struct.TransmitAntennaArray = MergeStructs(...
    ant_struct.TransmitAntennaArray,cdl_struct.TransmitAntennaArray);

ant_struct.ReceiveAntennaArray = MergeStructs(...
    ant_struct.ReceiveAntennaArray,cdl_struct.ReceiveAntennaArray);

cdl_struct.TransmitAntennaArray = ant_struct.TransmitAntennaArray;
cdl_struct.ReceiveAntennaArray  = ant_struct.ReceiveAntennaArray;

end

function s1 = MergeStructs(s1,s2)
%MERGESTRUCTS Merge sub-struct s2's fields and values into main struct s1
    f = fieldnames(s2);
    for i=1:length(f)
        s1.(f{i}) = s2.(f{i});
    end
end