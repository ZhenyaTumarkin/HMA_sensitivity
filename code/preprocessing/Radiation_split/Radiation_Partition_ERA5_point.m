%%%%%%%%%%%%%%%%%%%%%%
%Adapted from Achille's 'Radiatn_Partition_ERA5_arr'
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% close all; clc;delete(gcp('nocreate'));
tic

% glacier_id = 'RGI60-15.07886'; defined in bash script
data_root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing';
code_root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing';
point_path =  strcat(data_root,'/All_glaciers/',glacier_id,'/coords_out_',glacier_id,'.csv');
points = readtable(point_path);
num_points = height(points);

% glacier_id = 'RGI60-15.07886';

% x1=datetime('01-Oct-2019 00:00:00');
% x2=datetime('30-Sep-2020 23:00:00');



% ncor = feature('numcores');
% disp(ncor)
% parpool('local',4);


points_data = readtable(strcat(data_root, '/All_glaciers/',glacier_id, '/coords_out_',glacier_id,'.csv'));

points_length=size(points_data,1);
points_all = 1:num_points;




% Add path for meteorological forcing and DEM


        
forcing_location = strcat(data_root,'All_glaciers/',glacier_id);
addpath(genpath(strcat(data_root,forcing_location)))
%addpath(genpath([code_root 'Radiation_split/TC_Preprocessing_Functions']))
addpath(genpath([code_root 'Radiation_split']))




for point_num = points_all
%delta GMT from longitude (assume all points in same timezone - hence only
point_lon = (360-points.lon(point_num));   % DEGREES WEST!
DeltaGMT=timezone(point_lon);


Lon = points_data.lon(point_num);
Lat = points_data.lat(point_num);
elev = points_data.elev_m(point_num);


% disp(['Partitionning SWin at ' glacier_id ' at point ' string(point_num) ' from ' datestr(x1,'dd-mmm-yyyy HH:MM') ' to '...
%     datestr(x2,'dd-mmm-yyyy HH:MM')])

era5_data = parquetread(strcat(data_root,'/All_glaciers/',glacier_id,'/ERA5L_RH/', string(point_num), '.parquet'));


Date = era5_data.time;
daysPerMonth = eomday(year(Date),month(Date));

[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO; Datam(:,5)= MI; Datam(:,6)= SE;



PP = era5_data.tp;
Tdew = era5_data.d2m;
Swin = era5_data.ssrd;

% Use downscaled or bias-corrected variable
t=1;
N_time_step=length(Datam);
Datam_start = Datam(t,:);

YR = Datam_start(1);
MO = Datam_start(2);

nHours = size(Swin,1);
Datam_end = datetime(Datam_start) + hours(nHours-1);
DateHR = datetime(Datam_start):hours(1):datetime(Datam_end);


[SAD1,SAD2,SAB1,SAB2,PARB,PARD]=Automatic_Radiation_Partition_fast(datenum(DateHR),Lat,Lon,elev,DeltaGMT,squeeze(PP),squeeze(Tdew)-273.15,squeeze(Swin),1,0);
SAD1=SAD1.';
SAD2 = SAD2.';
SAB1 = SAB1.';
SAB2 = SAB2.';
PARB = PARB.';
PARD = PARD.';

LWIN = era5_data.strd;

TA  = era5_data.t2m  - 273.15;
Sp = era5_data.sp / 100;
Ws = era5_data.Ws;
time = era5_data.time;
RH = era5_data.RH;


out_data =  table(time,PP,Ws,Sp,LWIN,RH,TA,SAD1,SAD2,SAB1,SAB2,PARB,PARD);
parquetwrite(strcat(data_root,'/All_glaciers/',glacier_id,'/ERA5L_radiation_partitioned/',string(point_num), '.parquet'),out_data)
%writetable(out_data,strcat(data_root,'/All_glaciers/',glacier_id,'/ERA5L_radiation_partitioned/',string(point_num), '.dat'),'WriteRowNames',true)
end
% T = table(string(DateHR(:)), 'VariableNames', {'Date_time'});
% writetable(T, strcat(root, 'All_glaciers/', glacier_id, '/ERA5L_radiation_partitioned/time.csv'));


toc
