
clear;clc;close all

% --------------- changelog -------------
% v0.01 added sample capture switch
% v0.02 added 3rd degree polynomial fitting output of adjacent lanepoints 

% % _x, _y  = 일반차선
% % _cx,_cy = 중앙선
% % _bx,_by = 경계선

% % tar_ : 위경도 좌표계 / 회전 미적용 
% %       / 자차기준 단순 roi로 추출 (가로 80m,세로60m) / +y : 정북쪽

% % tar_ _bf : m 좌표계 / 자차기준 회전 적용
% %           / +y : 좌측 횡방향

% % tar2_ _bf : m좌표계 / 자차기준 회전 적용
% %             / tar__bf 에 차선번호 계산용 roi 적용 (종방향 5m / 횡방향 60m)

% % tar3_ _bf : m좌표계 / 자차기준 회전 적용
% %             / tar2__bf 중에서 현재 진행방향 차선 점만 추출 (중앙선 우측) - c,b 존재하지않음



% =============== TYPE FILE PATH AND NAME ===================

kml_file_path='C:\Users\jsss1\Desktop\ACL\pangyo_lane_annotation\kodas_lane_annotation_code\input'; % folder/folder/folder
kml_file_name='pangyo_a1_lane_kml.kml'; % '<filename>.kml'

kodas_data_path='C:\Users\jsss1\Desktop\ACL\pangyo_lane_annotation\kodas_lane_annotation_code\input'; % folder/folder/folder
kodas_data_name='kodas_example.mat';  % '<filename>.mat'
switch_save_movie=0

switch_plot_certain_sample=0
target_sample_index=2618;



% ================= Plot parameter==========================

roi_x1=100;
roi_y1=100; % too small value can cause lane num calc error
% plot 범위 설정 - viewrange 검색


% ============================================================





kmlpath=[kml_file_path '/' kml_file_name]
kodaspath=[kodas_data_path '/' kodas_data_name]

% ------------ Preprocessing (run once per kml) ----------------
% kml 2 mat
cd('./1.kml2struct')
ks18=kml2struct_surface_18(kmlpath);
save('../2.save_surface_coord/ks18mat.mat');
clearvars -except kodaspath switch* target* roi*
cd('../')

% extract lane points
cd('./2.save_surface_coord')
save_surface_byidx_18('ks18mat.mat','lanepoints18.mat')
cd('../')




% -------------------- lane num annotation -----------------------
clearvars -except kodaspath switch* target* roi*


load(kodaspath);
load('./3.lane_annotation/lanepoints18.mat');
total_lane_array=[];
cur_lane_array=[];
left_offset_array=[];
right_offset_array=[];
l_coeff_array=[];
r_coeff_array=[];


            if switch_save_movie
            close all

            AVi_File_Name = ['lanefitting_test' '.avi' ];
            disp('RECORDING...');
            aviobj = VideoWriter(AVi_File_Name);
            hFig=figure();
            end


for i=1:1:length(longitude)
    if switch_plot_certain_sample && i~=target_sample_index
        continue
    end
    cur_x=longitude(i);
    cur_y=latitude(i);
    cur_yaw=pi/2-IMU_Yaw(i); % KODAS data has different polar coordinates axis ( NORTH:0 // EAST:+pi/2 )
                             % This code uses ( NORTH : +pi/2 // EAST : 0 )
                             


VHCL_WIDTH=1.805;
VHCL_LENGTH=4.375;

l_mind_cand=[];
r_mind_cand=[];


% % tar_ : 위경도 좌표계 / 회전 미적용 
% %       / 자차기준 단순 roi로 추출 (가로 80m,세로60m) / +y : 정북쪽
        % normal line (white)
tar_x=total_long((-roi_x1/2/111320)<total_long-cur_x&total_long-cur_x<(+roi_x1/2/111320)&(-roi_y1/2/111320)<total_lat-cur_y&total_lat-cur_y<(+roi_y1/2/111320));
tar_y=total_lat((-roi_x1/2/111320)<total_long-cur_x&total_long-cur_x<(+roi_x1/2/111320)&(-roi_y1/2/111320)<total_lat-cur_y&total_lat-cur_y<(+roi_y1/2/111320));
        % center line              
tar_cx=total_long_c((-roi_x1/2/111320)<total_long_c-cur_x&total_long_c-cur_x<(+roi_x1/2/111320)&(-roi_y1/2/111320)<total_lat_c-cur_y&total_lat_c-cur_y<(+roi_y1/2/111320));
tar_cy=total_lat_c((-roi_x1/2/111320)<total_long_c-cur_x&total_long_c-cur_x<(+roi_x1/2/111320)&(-roi_y1/2/111320)<total_lat_c-cur_y&total_lat_c-cur_y<(+roi_y1/2/111320));
        % boundary line (yellow)
tar_bx=total_long_b((-roi_x1/2/111320)<total_long_b-cur_x&total_long_b-cur_x<(+roi_x1/2/111320)&(-roi_y1/2/111320)<total_lat_b-cur_y&total_lat_b-cur_y<(+roi_y1/2/111320));
tar_by=total_lat_b((-roi_x1/2/111320)<total_long_b-cur_x&total_long_b-cur_x<(+roi_x1/2/111320)&(-roi_y1/2/111320)<total_lat_b-cur_y&total_lat_b-cur_y<(+roi_y1/2/111320));



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

% convert to m scale  
tar_x_bf=tar_x_bf*111320;
tar_y_bf=tar_y_bf*111320;
tar_cx_bf=tar_cx_bf*111320;
tar_cy_bf=tar_cy_bf*111320;
tar_bx_bf=tar_bx_bf*111320;
tar_by_bf=tar_by_bf*111320;

tar_xtotal_bf=[tar_cx_bf;tar_x_bf;tar_bx_bf];
tar_ytotal_bf=[tar_cy_bf;tar_y_bf;tar_by_bf];

% car box on bodi-fixed coord.
bxh=[-VHCL_LENGTH+1.7774 1.7774 1.7774 -VHCL_LENGTH+1.7774]; byh=[-VHCL_WIDTH/2 -VHCL_WIDTH/2 VHCL_WIDTH/2 VHCL_WIDTH/2];
bbox_vert=polyshape(byh,bxh);

% 2nd ROI filtering to calculate cur_lane/total_lane

    % normal line (white)
tar2_x_bf=tar_x_bf((-2.5)<tar_x_bf&tar_x_bf<(+2.5)&(-30)<tar_y_bf&tar_y_bf<(+30));
tar2_y_bf=tar_y_bf((-2.5)<tar_x_bf&tar_x_bf<(+2.5)&(-30)<tar_y_bf&tar_y_bf<(+30));
if ~isempty(tar2_y_bf)
if sum(tar2_y_bf>0)>0
l_mind_cand(end+1,1)=min(tar2_y_bf(tar2_y_bf>0));
end
if sum(tar2_y_bf<0)>0
r_mind_cand(end+1,1)=max(tar2_y_bf(tar2_y_bf<0));
end
end
    % center line
tar2_cx_bf=tar_cx_bf((-2.5)<tar_cx_bf&tar_cx_bf<(+2.5)&(-30)<tar_cy_bf&tar_cy_bf<(+30));
tar2_cy_bf=tar_cy_bf((-2.5)<tar_cx_bf&tar_cx_bf<(+2.5)&(-30)<tar_cy_bf&tar_cy_bf<(+30));
if ~isempty(tar2_cy_bf)
if sum(tar2_cy_bf>0)>0
l_mind_cand(end+1,1)=min(tar2_cy_bf(tar2_cy_bf>0));
end
if sum(tar2_cy_bf<0)>0
r_mind_cand(end+1,1)=max(tar2_cy_bf(tar2_cy_bf<0));
end
end
    % boundary line (yellow)
tar2_bx_bf=tar_bx_bf((-2.5)<tar_bx_bf&tar_bx_bf<(+2.5)&(-30)<tar_by_bf&tar_by_bf<(+30));
tar2_by_bf=tar_by_bf((-2.5)<tar_bx_bf&tar_bx_bf<(+2.5)&(-30)<tar_by_bf&tar_by_bf<(+30));
if ~isempty(tar2_by_bf)
if sum(tar2_by_bf>0)>0
l_mind_cand(end+1,1)=min(tar2_by_bf(tar2_by_bf>0));
end
if sum(tar2_by_bf<0)>0
r_mind_cand(end+1,1)=max(tar2_by_bf(tar2_by_bf<0));
end
end


if ~isempty(l_mind_cand)
    l_offset=min(l_mind_cand);
else
    l_offset=NaN;
end
if ~isempty(r_mind_cand)
    r_offset=max(r_mind_cand);
else
    r_offset=NaN;
end

% find y value(lateral coord) on car position:x=0
if length(tar2_cy_bf)>0 && sum(tar2_cy_bf<0)==0
thc=min(tar2_cy_bf);

% extract normal line lanepoints on ongoing direction only
tar3_x_bf=tar2_x_bf(thc>tar2_y_bf);
tar3_y_bf=tar2_y_bf(thc>tar2_y_bf);
no_centerline_flag=0;
else
tar3_x_bf=[];
tar3_y_bf=[];
total_lane=NaN;
cur_lane=NaN;
no_centerline_flag=1;
end

if length(tar3_x_bf)>0
s=sort(tar3_y_bf);
ds=s-[s(1); s(1:end-1)];
sp=s(s>0);
total_lane=sum(ds>1.5)+2;
    if length(sp)>0
    dsp=sp-[sp(1); sp(1:end-1)];
    cur_lane=sum(dsp>1.5)+2;
    else
    cur_lane=1;
    end

elseif no_centerline_flag==1
total_lane=NaN;
cur_lane=NaN;
l_offset=NaN;
r_offset=NaN;

else 
total_lane=1;
cur_lane=1;
end





% grouping adjacent lanes

tar2_xtotal_bf=[tar2_cx_bf;tar2_x_bf;tar2_bx_bf];
tar2_ytotal_bf=[tar2_cy_bf;tar2_y_bf;tar2_by_bf];
l_fit_x=[]; l_fit_y=[];
r_fit_x=[]; r_fit_y=[];
l_lane_flag=0; r_lane_flag=0;
l_coeff=[0 0 0 0];
r_coeff=[0 0 0 0];

if ~isempty(tar2_xtotal_bf)
    
if sum((tar2_ytotal_bf<3)&(tar2_ytotal_bf>0.8))>0
    if sum((min(tar2_ytotal_bf(tar2_ytotal_bf>0.8))+0.3>tar2_ytotal_bf))>0
l_fit_x=tar2_xtotal_bf((tar2_ytotal_bf<3)&(tar2_ytotal_bf>0.8)&(min(tar2_ytotal_bf(tar2_ytotal_bf>0.8))+0.3>tar2_ytotal_bf));
l_fit_y=tar2_ytotal_bf((tar2_ytotal_bf<3)&(tar2_ytotal_bf>0.8)&(min(tar2_ytotal_bf(tar2_ytotal_bf>0.8))+0.3>tar2_ytotal_bf));
l_lane_flag=1;
    end
end
if sum((tar2_ytotal_bf>-3)&(tar2_ytotal_bf<-0.8))>0
    if sum((max(tar2_ytotal_bf(tar2_ytotal_bf<0.8))-0.3<tar2_ytotal_bf))>0
r_fit_x=tar2_xtotal_bf((tar2_ytotal_bf>-3)&(tar2_ytotal_bf<-0.8)&(max(tar2_ytotal_bf(tar2_ytotal_bf<0.8))-0.3<tar2_ytotal_bf));
r_fit_y=tar2_ytotal_bf((tar2_ytotal_bf>-3)&(tar2_ytotal_bf<-0.8)&(max(tar2_ytotal_bf(tar2_ytotal_bf<0.8))-0.3<tar2_ytotal_bf));
r_lane_flag=1;
    end
end
[~,cur_lfit_idx]=max(l_fit_x);
[~,cur_rfit_idx]=max(r_fit_x);

if isempty(l_fit_x)
    l_lane_flag=0;
end
if isempty(r_fit_x)
    r_lane_flag=0;
end

    % left-side 
    if l_lane_flag==1
    for j=1:20
    [~,cur_lfit_idx]=max(l_fit_x);
    tmpl=sqrt((tar_ytotal_bf-l_fit_y(cur_lfit_idx)).^2+(tar_xtotal_bf-l_fit_x(cur_lfit_idx)).^2/25);
    tmpl(tmpl==0)=100;
    tmpl(tar_xtotal_bf<=l_fit_x(cur_lfit_idx))=100;
    tmpl(abs(tar_ytotal_bf-l_fit_y(cur_lfit_idx))>1)=100;
    [~,next_lfit_idx]=min(tmpl);
    if min(tmpl)>=100
        break
    end
    l_fit_x=[l_fit_x; tar_xtotal_bf(next_lfit_idx)];
    l_fit_y=[l_fit_y; tar_ytotal_bf(next_lfit_idx)];
    end
    l_coeff=polyfit(l_fit_x,l_fit_y,3);
    end

    % right-side
    if r_lane_flag==1
    for k=1:20
    [~,cur_rfit_idx]=max(r_fit_x);
    tmpr=sqrt((tar_ytotal_bf-r_fit_x(cur_rfit_idx)).^2+(tar_xtotal_bf-r_fit_x(cur_rfit_idx)).^2/25);
    tmpr(tmpr==0)=100;
    tmpr(tar_xtotal_bf<=r_fit_x(cur_rfit_idx))=100;
    tmpr(abs(tar_ytotal_bf-r_fit_y(cur_rfit_idx))>1)=100;
    [~,next_rfit_idx]=min(tmpr);
    if min(tmpr)>=100
        break
    end
    r_fit_x=[r_fit_x; tar_xtotal_bf(next_rfit_idx)];
    r_fit_y=[r_fit_y; tar_ytotal_bf(next_rfit_idx)];
    end
    r_coeff=polyfit(r_fit_x,r_fit_y,3);
    end

end





% visualization
            if switch_save_movie || switch_plot_certain_sample
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
            plot(l_fit_y,l_fit_x,'og')
            plot(r_fit_y,r_fit_x,'om')
            tc=[ num2str(cur_lane) ' / ' num2str(total_lane) ];
            ty=['L:' string(l_coeff) 'R:' string(r_coeff) ];
            title(tc)
            yl=ylabel(ty);
            
            set(yl,'rotation',0,'VerticalAlignment','middle','FontSize',10);
            yl.Position(1)=55;
            yl.Position(2)=0;
            
% % viewrange
            ylim([-20 60])
            xlim([-40 40])
            grid on


            hold off
            

            end
            if switch_save_movie
            % video generation
            mo = getframe(hFig);
                    open(aviobj)
                    writeVideo(aviobj,mo);
                      drawnow

            end
  
    
total_lane_array(end+1,1)=total_lane;
cur_lane_array(end+1,1)=cur_lane;
left_offset_array(end+1,1)=l_offset;
right_offset_array(end+1,1)=r_offset;
l_coeff_array(end+1,:)=l_coeff;
r_coeff_array(end+1,:)=r_coeff;
end

            if switch_save_movie
            disp('Video Save End');
                close(aviobj);
                close(hFig);
            end

if ~switch_plot_certain_sample
clearvars tar* thc s sp VHCL* cur_lane total_lane cur_x cur_y cur_yaw bbox_vert bxh byh d* i total_lat* total_long* *_mind_cand *_offset roi_* switch* tmp* *_coeff *_flag *_fit_x *_fit_y j k *_idx
figure
plot(total_lane_array,'LineWidth',2)
hold on
plot(cur_lane_array,'LineWidth',2)
legend('total lane','cur lane')
ylim([0 7])
end
