function phase_struct = GenerateInitialPhase(cdl)
switch lower(cdl.RandomStream)
    case 'mt19937ar with seed'
        randStream = RandStream('mt19937ar',Seed=cdl.Seed);
    otherwise
        rng(cdl.Seed)
end

c0 = 299792458; % Light speed

phase_struct.InitPhase = cdl.InitialPhases;


