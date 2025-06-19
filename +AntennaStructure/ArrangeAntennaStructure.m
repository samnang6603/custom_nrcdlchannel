function mdl_out = ArrangeAntennaStructure( ...
    TransmitAntennaArray,...
    ReceiveAntennaArray,...
    TransmitArrayOrientation,...
    ReceiveArrayOrientation,...
    CarrierFrequency)

mdl_out.NumInputSignals = prod(TransmitAntennaArray.Size);
mdl_out.NumOutputSignals = prod(ReceiveAntennaArray.Size);

% Arrange Tx antenna layout
txLayout = AntennaStructure.ArrangeLayout(TransmitAntennaArray,TransmitArrayOrientation);
txSize = size(txLayout.ElementPositions,2:6); % idx 2-6 to take into account multiple panels
txSubarraySize = size(txLayout.SubarrayPositions);
txElemPattern = AntennaStructure.ElementPattern.GenerateElementPattern();
mdl_out.TransmitAntennaArray = AllocateArrayStruct(txSize,txSubarraySize,txElemPattern);

% Arrange Rx antenna layout
rxLayout = AntennaStructure.ArrangeLayout(ReceiveAntennaArray,ReceiveArrayOrientation);
rxSize = size(rxLayout.ElementPositions,2:6); % idx 2-6 to take into account multiple panels
rxSubarraySize = size(rxLayout.SubarrayPositions);
rxElemPattern = AntennaStructure.ElementPattern.GenerateElementPattern();
mdl_out.ReceiveAntennaArray = AllocateArrayStruct(rxSize,rxSubarraySize,rxElemPattern);

% Merge Tx/Rx layout to main struct for output
mdl_out.TransmitAntennaArray = MergeStructs(mdl_out.TransmitAntennaArray,txLayout);
mdl_out.ReceiveAntennaArray = MergeStructs(mdl_out.ReceiveAntennaArray,rxLayout);

end

function mdl_struct = AllocateArrayStruct(asize,subarraySize,elemPattern)
% ALLOCATEARRAYSTRUCT Allocate empty sub-struct to be merged with main 
% struct
mdl_struct = struct('Position',zeros([3 1]),...
    'InitialArrayOrientation',zeros(3,1),...
    'Orientation',zeros([3 1]),...
    'ElementPositions',zeros([3 asize]),...
    'SubarrayPositions',zeros(subarraySize),...
    'ElementOrientations',zeros([3 asize]),...
    'ElementPattern',elemPattern);
end

function s1 = MergeStructs(s1,s2)
%MERGESTRUCTS Merge sub-struct s2's fields and values into main struct s1
    f = fieldnames(s2);
    for i=1:length(f)
        s1.(f{i}) = s2.(f{i});
    end
end