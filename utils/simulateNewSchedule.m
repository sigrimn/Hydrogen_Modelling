function wellsol = simulateNewSchedule(scheduleRef, wellRef, setup, dt)

% schedule = simpleSchedule(dt, 'W', well);
% 
% %Setting controls:
% for kw = 1:scheduleRef.step.control(end)
%     schedule.control(kw).W.val = scheduleRef.control(kw).W.val;
%     %schedule.control(kw).bc = scheduleRef.control(1).bc;
% end
% 
% wellsol = simulateScheduleAD(setup.state0, setup.model, schedule);


wOpt = setup.schedule.control(1).W;

schedule = simpleSchedule(dt, 'W', wOpt);

%Setting controls:
for kw = 1:numel(wellRef)
    schedule.control.W(kw).val = scheduleRef.control.W(kw).val;
end

wellsol = simulateScheduleAD(setup.state0, setup.model, schedule);
