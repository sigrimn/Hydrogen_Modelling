function [scheduleCoarse, stateCoarse0] = makeCoarseSchedule3D(modelCoarse, modelRef, scheduleRef, state0)

stateCoarse0   = upscaleState(modelCoarse, modelRef, state0);

scheduleCoarse = upscaleSchedule(modelCoarse, scheduleRef);

%Validating the schedule:
scheduleCoarse = validateSchedule(modelCoarse, scheduleCoarse);
