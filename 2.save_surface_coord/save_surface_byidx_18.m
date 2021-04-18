function save_surface_byidx_18(kmlmatname,lanepointsmatname)


% =============================
%  I : pangyo_ks18.mat
%  O : pangyo_ks_var_18.mat
% =============================

% load('pangyo_kml_2_mat.mat')
load(kmlmatname)
% lat_init = 37.58;
% long_lnit = 126.89;




total_lane_rlinkid=[];total_lane_llinkid=[];total_lane_id=[];
total_lane_lanecode=[];total_lane_lanetype=[];total_lane_lno=[];
total_lane_flag_c=zeros(length(ks18),1);
total_long=[];total_lat=[];
total_long_c=[];total_lat_c=[];
total_long_b=[];total_lat_b=[];

tmprx=[];tmpry=[];tmprx2=[];tmpry2=[];tmprx3=[];tmpry3=[];tmprx4=[];tmpry4=[];tmprx5=[];tmpry5=[];
count=0;count2=0;count3=0;count4=0;count5=0;
% for i=1:1:length(section)
% n=section(i);

for j=1:1:length(ks18)
    % n=1;
    % for j=1:20
    cur_lat=ks18(j).Lat;
    cur_long=ks18(j).Lon;   
    cur_rlinkid=ks18(j).r_linkid;
    cur_llinkid=ks18(j).l_linkid;
    cur_lanecode=ks18(j).lanecode;
    cur_id=ks18(j).id;
    cur_lno=ks18(j).lno;
    cur_lanetype=ks18(j).lanetype;
    tx=[];ty=[];
    if cur_lanecode=="01" || cur_lanecode=="02"
        total_lane_flag_c(j)=1;
        total_long_c=[total_long_c;cur_long];
        total_lat_c=[total_lat_c;cur_lat];
    elseif cur_lanecode=="03" || cur_lanecode=="05"
        total_long=[total_long;cur_long];
        total_lat=[total_lat;cur_lat];
    elseif cur_lanecode=="07" || cur_lanecode=="08" || cur_lanecode=="09"
        total_long_b=[total_long_b;cur_long];
        total_lat_b=[total_lat_b;cur_lat];
    end
    
 
    total_lane_rlinkid=[total_lane_rlinkid;cur_rlinkid];
    total_lane_llinkid=[total_lane_llinkid;cur_llinkid];
    total_lane_id=[total_lane_id;cur_id];
    total_lane_lanecode=[total_lane_lanecode;cur_lanecode];
    
   
    
    
%     
% 
%     for i=1:1:length(cur_lat)
%         [tx(i) ty(i)]=GaussKruegerProjection(cur_lat(i),cur_long(i),lat_init,long_lnit);
%     end
%     theta = 0.0085+7.203134235761466e-04+4.993072981579472e-04;
% 
%     tans=[];
%     tans2=[];
%     for i =1 : length(tx)
%         Trans_matrix = [cos(theta) sin(theta); -sin(theta) cos(theta)];
%         tans(i,1:2) = Trans_matrix * [tx(1,i) ty(1,i)]';
        tans2=[cur_long, cur_lat];
%     end
%     total_lx=[total_lx;tans2(:,1)];
%     total_ly=[total_ly;tans2(:,2)]; % long lat 2 xy coordinate done
    
    eval(['laneraw_' num2str(j) '=tans2;']) % save points group as ks18 index
    
% GEN MATRIX WITH NAME OF ITS ID
  
    
    


end
figure
plot(total_long,total_lat,'.k');
grid on;axis equal;hold on;
plot(total_long_c,total_lat_c,'.r')
plot(total_long_b,total_lat_b,'.b')

fullpath=['../3.lane_annotation/' lanepointsmatname];
save(fullpath,'total_lat*','total_long*');
clearvars -except total_lat* total_long* lanepointsmatname
