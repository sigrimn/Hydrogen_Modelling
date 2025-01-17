function [scheduleCoarse, stateCoarse0] = makeCoarseSchedule2D(modelCoarse, modelRef, scheduleRef, state0)

stateCoarse0   = upscaleState(modelCoarse, modelRef, state0);
% Changing cstatus which is which of the well cells that are active or not

stateCoarse0.wellSol.cstatus = logical([1, zeros(1, 7)]);

scheduleCoarse = upscaleSchedule(modelCoarse, scheduleRef);

%Validating the schedule:
scheduleCoarse = validateSchedule(modelCoarse, scheduleCoarse);