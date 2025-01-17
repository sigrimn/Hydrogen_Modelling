function [bhp_error, h2o_rate_error, h2_rate_error] = calculate_well_errors(ref_well_sols, tuned_well_sols)
    % Function to calculate errors between reference and tuned well solutions
    %
    % Inputs:
    %   ref_well_sols   - Cell array or array containing reference well solutions for all timesteps
    %   tuned_well_sols - Cell array or array containing tuned well solutions for all timesteps
    %
    % Outputs:
    %   bhp_error       - Total error in bottom-hole pressure (BHP) over all timesteps
    %   h2o_rate_error  - Total error in water rate (H2O rate) over all timesteps
    %   h2_rate_error   - Total error in hydrogen rate (H2 rate) over all timesteps

    % Initialize errors
    bhp_error = 0;
    h2o_rate_error = 0;
    h2_rate_error = 0;
    
    % Iterate through each timestep
    for t = 1:length(ref_well_sols)
        % Extract BHP, H2O rate, and H2 rate for the current timestep
        ref_bhp = ref_well_sols{t}.bhp;
        tuned_bhp = tuned_well_sols{t}.bhp;

        ref_h2o_rate = ref_well_sols{t}.qOs; % Water production rate
        tuned_h2o_rate = tuned_well_sols{t}.qOs;

        ref_h2_rate = ref_well_sols{t}.qGs; % Gas (hydrogen) production rate
        tuned_h2_rate = tuned_well_sols{t}.qGs;

        % Accumulate absolute errors for this timestep
        bhp_error = bhp_error + abs(ref_bhp - tuned_bhp);
        h2o_rate_error = h2o_rate_error + abs(ref_h2o_rate - tuned_h2o_rate);
        h2_rate_error = h2_rate_error + abs(ref_h2_rate - tuned_h2_rate);
    end
    bhp_error = bhp_error/length(ref_well_sols);
    h2o_rate_error = h2o_rate_error /length(ref_well_sols);
    h2_rate_error = h2_rate_error/length(ref_well_sols);
end
