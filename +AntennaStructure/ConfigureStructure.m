function ant_struct = ConfigureStructure(cdl_struct)

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
end