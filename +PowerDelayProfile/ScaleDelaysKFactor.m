function pdp_table_out = ScaleDelaysKFactor(pdp_table,K_desired,targetDelaySpread)
% SCALEDELAYSKFACTOR
%
% Initial set of delay profile idx
% TR 38.901 7.7-1 CDL PDP table for info on indices
[delayIdx,powIdx,AoD,AoA,ZoD,ZoA] = deal(1,2,3,4,5,6);

pdp_table_out = pdp_table; % save pdp_table for debugging

% If LoS exist (CDL-D, CDL-E), apply K-Factor scaling
if ~isnan(K_desired)

    % calculate K-Factor from pdp, lookup power index
    % See Equation 7.7.6-2
    total_pow_lin = sum(10.^(pdp_table_out(2:end,powIdx)/10)); % total power excluding LOS (0-delay) cluster
    K_model = pdp_table_out(1,powIdx) - pow2db(total_pow_lin); % subtract from LOS cluster

    % scale the power of all taps excluding the 1st tap (0-delay
    % component) by subtract the difference of target K-Factor and the
    % calculated K-Factor
    % See Equation 7.7.6-1
    P_model = pdp_table_out(2:end,powIdx);
    pdp_table_out(2:end,powIdx) =  P_model - K_desired + K_model;

    % calculate the RMS delay spread after the K-factor adjustment
    ptmp = 10.^(pdp_table_out(:,powIdx)/10);
    tau = pdp_table_out(:,delayIdx);
    w1 = sum(ptmp.*tau)/sum(ptmp); % weighted rms w1
    w2 = sum(ptmp.*tau.^2)/sum(ptmp); % weighted rms w2
    tau_rms = sqrt(w2 - w1^2);

    % normalize delay spread
    pdp_table_out(:,delayIdx) = pdp_table_out(:,delayIdx)/tau_rms;
end

% % scale delays according to desired delay spread, Section 7.7.3
pdp_table_out(:,delayIdx) = pdp_table_out(:,delayIdx)*targetDelaySpread;

end



