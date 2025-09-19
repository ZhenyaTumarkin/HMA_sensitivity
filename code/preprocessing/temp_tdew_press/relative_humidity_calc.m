%relative humidity
tic


root = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing';

addpath(genpath('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/preprocessing/temp_tdew_press/ConvertHumidity'))
%now load data
% Id = 'RGI60-15.07886';

point_path =  strcat(root,'/All_glaciers/',Id,'/coords_out_',Id,'.csv');
era5L_path = strcat(root,'/All_glaciers/',Id,'/ERA5L_at_gridpoint');

points = readtable(point_path);
num_points = height(points);
disp('number of points = '); disp(num_points);

for i = 1:num_points
    point_data = parquetread(strcat(era5L_path,'/',num2str(i),'.parquet'));
    point_data.RH  = convert_humidity (point_data.sp, point_data.t2m,  point_data.d2m, 'dew point temperature', 'relative humidity', 'Murphy&Koop2005', 0);
    parquetwrite(strcat(root,'/All_glaciers/',Id,'/ERA5L_RH/',num2str(i),'.parquet'),point_data)   % parquet is a much more compressed file type, about 5x faster to write 
    % writetable(point_data,strcat(root,'/All_glaciers/',Id,'/ERA5L_RH/',num2str(i),'.dat'))
end

disp('RH complete ')

toc

