f=cursor_info.Position

for i=1:1:length(total_lane_rlinkid)
a=eval(['sum(sum(laneraw_' num2str(i) '==f(1)));']);
if a>0
i
break
end
end
total_lane_lanecode(i)
