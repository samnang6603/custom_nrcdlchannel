function antstruct = ConfigureStructure(cdl)

TransmitAntennaArray = cdl.TransmitAntennaArray;
ReceiveAntennaArray = cdl.ReceiveAntennaArray;
TransmitArrayOrientation = cdl.TransmitArrayOrientation;
ReceiveArrayOrientation = cdl.ReceiveArrayOrientation;
CarrierFrequency = cdl.CarrierFrequency;

mdlArrays = AntennaStructure.ArrangeAntennaStructure(...
    TransmitAntennaArray,...
    ReceiveAntennaArray,...
    TransmitArrayOrientation,...
    ReceiveArrayOrientation,...
    CarrierFrequency);

antstruct = mdlArrays;
end