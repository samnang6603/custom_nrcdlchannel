function cdl_struct = AdvanceTime(cdl_struct)
cdl_struct.CurrentTime = cdl_struct.CurrentTime + ...
    cdl_struct.NumTimesSamples/cdl_struct.SampleRate;
end