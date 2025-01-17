function weighting = objectiveWeighting(wellSolsFine)
    % Calculate weights for the mismatch function 
    % based on fine-scale solution.
    
    qGs = getWellOutput(wellSolsFine, 'qGs');
    qOs = getWellOutput(wellSolsFine, 'qOs');
    bhp = getWellOutput(wellSolsFine, 'bhp');
    
    wqGs = 1/max(abs(qGs(:)));
    wqOs = 1/max(abs(qOs(:)));
    wbhp = 1/(max(bhp(:)) - min(bhp(:)));
    
    weighting =  {'WaterRateWeight',  wqGs, ...
                  'OilRateWeight',    wqOs, ...
                  'BHPWeight',        wbhp};
end
