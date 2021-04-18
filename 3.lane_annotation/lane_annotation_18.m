function [cur_lane,total_lane]=lane_annotation_18(cur_x,cur_y,cur_yaw,kml)
load(kml);




% load('CTRL_210128_006')
close all

%  127.1054  37.4024

% load daq data ( .mat format converted from .bag )
% cx=class_pos_gt.CarTx;
% cy=class_pos_gt.CarTy;  
% cyaw=class_vehicleinfo.VehicleYaw; 
VHCL_WIDTH=1.805/90000;
VHCL_LENGTH=4.375/90000;
% 
% % 
% % AVi_File_Name = ['laneinfo_recog' '.avi' ];
% % disp('RECORDING...');//
% % aviobj = VideoWriter(AVi_File_Name);
% % hFig=figure();
% 
% % set current timestamp (20hz)
% t=120;
% idx=20*t;
% % for idx=1:1:length(cx)
% % for idx=1161
% cur_x=cx(idx);  cur_y=cy(idx);  %cur_yaw=cyaw(idx);
% cur_x=127.1054; cur_y=37.4024; cur_yaw=116*pi/180;


% tunnel
% cur_x=127.10438; cur_y=37.4045; cur_yaw=135*pi/180;


tar_x=total_long((-40/90000)<total_long-cur_x&total_long-cur_x<(+40/90000)&(-30/90000)<total_lat-cur_y&total_lat-cur_y<(+30/90000));
tar_y=total_lat((-40/90000)<total_long-cur_x&total_long-cur_x<(+40/90000)&(-30/90000)<total_lat-cur_y&total_lat-cur_y<(+30/90000));
tar_cx=total_long_c((-40/90000)<total_long_c-cur_x&total_long_c-cur_x<(+40/90000)&(-30/90000)<total_lat_c-cur_y&total_lat_c-cur_y<(+30/90000));
tar_cy=total_lat_c((-40/90000)<total_long_c-cur_x&total_long_c-cur_x<(+40/90000)&(-30/90000)<total_lat_c-cur_y&total_lat_c-cur_y<(+30/90000));
tar_bx=total_long_b((-40/90000)<total_long_b-cur_x&total_long_b-cur_x<(+40/90000)&(-30/90000)<total_lat_b-cur_y&total_lat_b-cur_y<(+30/90000));
tar_by=total_lat_b((-40/90000)<total_long_b-cur_x&total_long_b-cur_x<(+40/90000)&(-30/90000)<total_lat_b-cur_y&total_lat_b-cur_y<(+30/90000));
% 
% 
% tar_x=total_long;
% tar_y=total_lat;
% tar_cx=total_long_c;
% tar_cy=total_lat_c;
% tar_bx=total_long_b;
% tar_by=total_lat_b;



 % rotate to body-fixed coord.
dx=(tar_x-cur_x);dy=(tar_y-cur_y);
dxc=(tar_cx-cur_x);dyc=(tar_cy-cur_y);
dxb=(tar_bx-cur_x);dyb=(tar_by-cur_y);
tar_x_bf=dx*cos(-cur_yaw)-dy*sin(-cur_yaw);
tar_y_bf=dx*sin(-cur_yaw)+dy*cos(-cur_yaw);
tar_cx_bf=dxc*cos(-cur_yaw)-dyc*sin(-cur_yaw);
tar_cy_bf=dxc*sin(-cur_yaw)+dyc*cos(-cur_yaw);
tar_bx_bf=dxb*cos(-cur_yaw)-dyb*sin(-cur_yaw);
tar_by_bf=dxb*sin(-cur_yaw)+dyb*cos(-cur_yaw);


% car box on bodi-fixed coord.
bxh=[-VHCL_LENGTH+1.7774/90000 1.7774/90000 1.7774/90000 -VHCL_LENGTH+1.7774/90000]; byh=[-VHCL_WIDTH/2 -VHCL_WIDTH/2 VHCL_WIDTH/2 VHCL_WIDTH/2];
bbox_vert=polyshape(byh,bxh);

% 2nd ROI filtering to calculate cur_lane/total_lane
tar2_x_bf=tar_x_bf((-2.5/90000)<tar_x_bf&tar_x_bf<(+2.5/90000)&(-30/90000)<tar_y_bf&tar_y_bf<(+30/90000));
tar2_y_bf=tar_y_bf((-2.5/90000)<tar_x_bf&tar_x_bf<(+2.5/90000)&(-30/90000)<tar_y_bf&tar_y_bf<(+30/90000));

tar2_cx_bf=tar_cx_bf((-2.5/90000)<tar_cx_bf&tar_cx_bf<(+2.5/90000)&(-30/90000)<tar_cy_bf&tar_cy_bf<(+30/90000));
tar2_cy_bf=tar_cy_bf((-2.5/90000)<tar_cx_bf&tar_cx_bf<(+2.5/90000)&(-30/90000)<tar_cy_bf&tar_cy_bf<(+30/90000));

tar2_bx_bf=tar_bx_bf((-2.5/90000)<tar_bx_bf&tar_bx_bf<(+2.5/90000)&(-30/90000)<tar_by_bf&tar_by_bf<(+30/90000));
tar2_by_bf=tar_by_bf((-2.5/90000)<tar_bx_bf&tar_bx_bf<(+2.5/90000)&(-30/90000)<tar_by_bf&tar_by_bf<(+30/90000));
% 
% thc=min(tar2_cy_bf);
% if thc<0
%     cur_yaw=cur_yaw+3.14;       
%    
% 
% 
%          % rotate to body-fixed coord.
%         dx=(tar_x-cur_x);dy=(tar_y-cur_y);
%         dxc=(tar_cx-cur_x);dyc=(tar_cy-cur_y);
%         dxb=(tar_bx-cur_x);dyb=(tar_by-cur_y);
%         tar_x_bf=dx*cos(-cur_yaw)-dy*sin(-cur_yaw);
%         tar_y_bf=dx*sin(-cur_yaw)+dy*cos(-cur_yaw);
%         tar_cx_bf=dxc*cos(-cur_yaw)-dyc*sin(-cur_yaw);
%         tar_cy_bf=dxc*sin(-cur_yaw)+dyc*cos(-cur_yaw);
%         tar_bx_bf=dxb*cos(-cur_yaw)-dyb*sin(-cur_yaw);
%         tar_by_bf=dxb*sin(-cur_yaw)+dyb*cos(-cur_yaw);
% 
% 
%         % car box on bodi-fixed coord.
%         bxh=[-VHCL_LENGTH+1.7774 1.7774 1.7774 -VHCL_LENGTH+1.7774]; byh=[-VHCL_WIDTH/2 -VHCL_WIDTH/2 VHCL_WIDTH/2 VHCL_WIDTH/2];
%         bbox_vert=polyshape(byh,bxh);
% 
%         % 2nd ROI filtering to calculate cur_lane/total_lane
%         tar2_x_bf=tar_x_bf((-2.5)<tar_x_bf&tar_x_bf<(+2.5)&(-30)<tar_y_bf&tar_y_bf<(+30));
%         tar2_y_bf=tar_y_bf((-2.5)<tar_x_bf&tar_x_bf<(+2.5)&(-30)<tar_y_bf&tar_y_bf<(+30));
% 
%         tar2_cx_bf=tar_cx_bf((-2.5)<tar_cx_bf&tar_cx_bf<(+2.5)&(-30)<tar_cy_bf&tar_cy_bf<(+30));
%         tar2_cy_bf=tar_cy_bf((-2.5)<tar_cx_bf&tar_cx_bf<(+2.5)&(-30)<tar_cy_bf&tar_cy_bf<(+30));
% 
%     
% end

% extract nearest center line threshold
if length(tar2_cy_bf)>0
thc=min(tar2_cy_bf);
tar3_x_bf=tar2_x_bf(thc>tar2_y_bf);
tar3_y_bf=tar2_y_bf(thc>tar2_y_bf);
else
tar3_x_bf=[];
tar3_y_bf=[];
end

if length(tar3_x_bf)>0
s=sort(tar3_y_bf);
ds=s-[s(1); s(1:end-1)];
sp=s(s>0);
total_lane=sum(ds>1.5/90000)+2;
    if length(sp)>0
    dsp=sp-[sp(1); sp(1:end-1)];
    cur_lane=sum(dsp>1.5/90000)+2;
    else
    cur_lane=1;
    end

else
total_lane=NaN;
cur_lane=NaN;
end
cla
hold on
axis equal
set(gca, 'XDir', 'reverse')
plot(bbox_vert)
plot(tar_y_bf,tar_x_bf,'.k')
plot(tar_cy_bf,tar_cx_bf,'.b')
plot(tar_by_bf,tar_bx_bf,'.m')

plot(tar2_y_bf,tar2_x_bf,'ok')
plot(tar2_cy_bf,tar2_cx_bf,'ob')


plot(tar3_y_bf,tar3_x_bf,'.r')
tc=[ num2str(cur_lane) ' / ' num2str(total_lane)];
title(tc)
ylim([-20/90000 60/90000])
xlim([-40/90000 40/90000])
grid on


hold off
% 
% mo = getframe(hFig);
%         open(aviobj)
%         writeVideo(aviobj,mo);
%           drawnow

  

% disp('Video Save End');
%     close(aviobj);
%     close(hFig);

