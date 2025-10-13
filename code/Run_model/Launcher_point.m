% feature('DisableExternalDDUX', 1);   %%%%stop sending data to mathworks (can apparently stop some crashes)
%%%%%%%% LINES TO BE CHANGED %%%%%%%%

model_path = '/nfs/scistore18/pelligrp/etumarki/TeC_code'; % Put here the path of where you downloaded the repository
run_path   = '/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/code/Run_model';


% glacier_id = 'RGI60-13.19847';
% point_id = 9; 
%precip_tune = 0.5;   %factor to multiply precipitation by
point_id = str2num(point_id);
data_path  = strcat('/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/preprocessing/All_glaciers/',glacier_id);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read in precipitation modification
tp_calib_path = ['/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/Outputs/' glacier_id '/tp_calib.csv'];
if isfile(tp_calib_path)
    p_mod_file = readtable(tp_calib_path);
    disp(p_mod_file)
    precip_tune= p_mod_file.next_p(end);
else
    precip_tune=1;
end
disp(['precipitation factor ' num2str(precip_tune)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study_name = 'single_point_run';

% How to store outputs %%%%
output_daily = 1;  % Recommended for long-term simulations (> 10-20 years, otherwise data volume is too important)


%%%
glacier_points = readtable(strcat(data_path,'/coords_out_',glacier_id,'.csv'));

Lat = glacier_points.lat(point_id);
Lon = glacier_points.lon(point_id);

% Study site details
FORCING = "ERA5Land";


%%%%%choose utm zone

minx = floor(Lon);
rem = mod(minx,6);
utm_border = (minx-rem)/6;
UTM_zone = int32(utm_border+31);


%DeltaGMT= 5;   % automate later
%delta GMT from longitude (assume all points in same timezone - hence only
point_lon = (360-glacier_points.lon(point_id));   % DEGREES WEST!
DeltaGMT=timezone(point_lon);

t_bef = 0.5; t_aft = 0.5;
  
% Select simulation period (start and end)
% x1s =  "01-Jul-2021 00:00:00"; % Starting point of the simulation
% % x2s =  "30-Sep-2023 23:00:00"; % Last timestep of the simulation
% date_start = datetime(x1s);
% date_end = datetime(x2s);

elevation = glacier_points.elev_m(point_id);

% Tmod = 0; % temperature modification above clean ice [°C];
% Pmod = 0; % factor Pmod preciptation modification, e.g. 0.3 means 30% more precipitation at highest elevation

%%%%%should read in as metadata!!!!
% Z_min = 3370; % lowest elevation for linear precipitation modification (min factor -> 0)
% Z_max = 5000; % highest elevation for linear precipitation modification (max factor -> Pmod)
%%%%%%%%%%%5


%%% Precipitation phase parametrization

%1 = 2-threshold, 2 = Ding 2017, 3 = single-threshold, 4 = Pomeroy 2013, 5
%= Wang 2019, 6 = Jennings 2018

parameterize_phase.OPT_Pr_Part = 2; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = 2; % Upper air temperature for dual temperature threshold
parameterize_phase.Tmin = 0; % Lower air temperature for dual temperature threshold
parameterize_phase.Tconst = 2; % Air temperature for constant thresholds

% Skin layer thickness:
hSTL = 0.003; %m

% Albedo scheme choice
Albsno_method = 5; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017
 

% Create the directy where model outputs will be stored

if ~exist('outlocation','var')
    disp('standard_out_loc chosen: data/Outputs/')
    outlocation = ['/nfs/scistore18/pelligrp/etumarki/HMA_sensitivity/data/Outputs/'  glacier_id '/run_' num2str(precip_tune,'%.3f') ];
end
if ~exist(outlocation, 'dir'); mkdir(outlocation); addpath(genpath(outlocation)); end
out = strcat(outlocation,'/point_run_',num2str(point_id),'.mat');%file path initial conditions
disp(outlocation)

%dependencies
addpath(genpath([run_path,'/Inputs/'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([model_path, '/Inputs'])); % Add path to Ca_Data
addpath(genpath([model_path, '/T&C_Code'])); % Add path to T&C codes



% Load CO2 data
load('Ca_Data.mat');
Ca_all = Ca;
topo = 0;    % include in future?????

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Impose measured albedo on glacier areas %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fn_alb_elev = [SITE '_Albedo_vs_elev.mat'];

% if exist(fn_alb_elev,'file')>0
% disp('Using measured glacier albedo')
% load([SITE '_Albedo_vs_elev'])

% Afirn = DTM.*0 + 0.28;
% dem_inc = DTM <= Alb_el_tt{1,2};
% Afirn(dem_inc) = Alb_el_tt{1,1};

% for ii = 1:size(Alb_el_tt,1)
%    if ii == size(Alb_el_tt,1)
%        dem_inc = DTM >= Alb_el_tt{ii,2};
%        Afirn(dem_inc) = Alb_el_tt{ii,1};
%    else
%        dem_inc = DTM > Alb_el_tt{ii,2} & DTM < Alb_el_tt{ii+1,2};
%        Afirn(dem_inc) = 0.5*(Alb_el_tt{ii,1} + Alb_el_tt{ii+1,1});
%    end 
% end 
% Afirn(Afirn > 0.4) = 0.4; %Limit bare ice albedo to 0.4, as above it's firn.
% else 
%     Afirn = DTM.*0 + 0.28;   
% end

%POIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POI = readtable([SITE '_MultiPoints.txt']); %import table with points info
% [POI.LAT, POI.LON] = utm2ll(POI.LON_UTM, POI.LAT_UTM, UTM_zone);

% % Select the indices of POI to run T&C for
% if exist('loc','var') %Option 1: if run on a cluster as an array task, use cluster array task ID as loc
%     LOC = loc; 
% else
%     LOC = 2; %1:size(POI,1); %loc = 11;       %Option 2: Manually select the indices
% end 

% for loc = LOC  
 
%get location info
% id_location = char(string(POI.Name(loc)));
% Lat = POI.LAT(loc);
% Lon = POI.LON(loc);




% Zbas = DTM(POI.ROW(loc),POI.COL(loc));
Zbas = glacier_points.elev_m(point_id);

%location name and cell index
% ij = POI.idx(loc);
% [i, j] = ind2sub(size(DTM), ij);


%% FORCING

forcing_all = parquetread(strcat(data_path,'/ERA5L_radiation_partitioned/',num2str(point_id),'.parquet'),VariableNames=["time","PP","Ws","Sp","Lwin","RH","Ta","SAD1","SAD2","SAB1","SAB2","PARB","PARD"]); % Load forcing table for the current POI
Date_all = forcing_all.time;
if ~exist('date_start','var')
    date_start = Date_all(1);
end
if ~exist('date_end','var')
    date_end = Date_all(end);
end




disp(['Site selected: ' glacier_id])
disp(['Forcing selected: ' char(FORCING)])
disp(['Running T&C for pixel: ' num2str(point_id)])
disp(['Simulation period: ' datestr(date_start) ' to ' datestr(date_end)])


% Time step 
dt=3600; %%[s]
dth=1; %%[h]
%define period and time zone info
x1=find(date_start == Date_all,1);
x2=find(date_end == Date_all,1);


Date = Date_all(x1:x2);
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

%load carbon data and narrow down to period
d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
% disp(d1);disp(d2);
Ca=Ca_all(d1:d2); 



clear d1 d2
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol]
%narrow down period of forcing data
forcing = forcing_all(x1:x2,:);
NN= height(forcing);%%% time Step

%height of virtual station
zatm_hourly = repmat(2.00,height(forcing),1);
zatm_surface = [18 18 2 2 18 18];
zatm_hourly_on=0;

%load all forcing data
Ameas = zeros(NN,1);

N=forcing.Lwin; Latm=forcing.Lwin;

% Precipitation
Pr=forcing.PP*precip_tune;Pr(isnan(Pr))=0; Pr(Pr<0.01)=0;

% if Pmod >0
%   Pr = Pr.*Pmod_S(ij);
% end


Pre=forcing.Sp;    
Ta=forcing.Ta; 

% Apply Tmod 
% if (GLH(ij)>0) && (DEB_MAP(ij) < 10) % glacier, but without debris
%     Ta(Ta >0) = Ta(Ta >0) + Tmod;
% end 

%%% Wind Speed
Ws=forcing.Ws*2.5 ; Ws(Ws < 0.05) = 0.05;
%%% Relative humidity
if max(forcing.RH)<= 1
    U=forcing.RH;
else
    U=forcing.RH./100;
end


%%% Radiation
SAD1=forcing.SAD1;SAD2=forcing.SAD2; SAB1=forcing.SAB1;SAB2=forcing.SAB2;
PARB=forcing.PARB; PARD=forcing.PARD;


alpha=0; %switch for albedo
Ameas_t=0; %albedo
Aice_meas_on_hourly=zeros(height(forcing),1); %albedo
Asno_meas_on_hourly=zeros(height(forcing),1); %albedo

%%% esat/ea/Ds/Tdew
a=17.27; b=237.3;
esat=611*exp(a*Ta./(b+Ta)); %Vapour pressure at saturation (Pa)
ea=U.*esat;                 %Vapour pressure (Pa)
Ds= esat - ea;              %Vapor Pressure Deficit (Pa)
Ds(Ds<0)=0; 
xr=a*Ta./(b+Ta)+log10(U);
Tdew=b*xr./(a-xr);               %Presumed dewpoint temperature (°C)
clear a b xr;

% Initial daily mean values for Ding parametrization

Ta_Ding_d = nanmean(Ta(1:24));
Pre_Ding_d = nanmean(Pre(1:24));
ea_Ding_d = nanmean(ea(1:24));

%% TOPOGRAPHY
% num_cell=numel(DTM);
% [m_cell,n_cell]=size(DTM);
% MASK = MASK.*0+1;
% MASKn=reshape(MASK,num_cell,1);

if topo == 1
    
    %load topography data and narrow down to period
    m = matfile([SITE '_ShF_' char(FORCING) '.mat']); % ShF matrix created during pre-processing step

    x1_top=find(date_start==m.Top_Date,1);
    x2_top=find(date_end==m.Top_Date,1); 

    ShF = double(squeeze(m.ShF_S(ij,:)));
    ShF = ShF(x1_top:x2_top)';
    
    rho_g = 0.35; %%% Spatial Albedo
    zeta_S = m.zeta_Sts(x1_top:x2_top,1);
    h_S = m.h_Sts(x1_top:x2_top,1);
    clear e1 e2 Top_Date zeta_Sts H_Sts ShF_S
    SvF = m.SvF_S(ij,1); % Sky view factor at pixel ij
    Ct = m.Ct_S(i,j);
    Slo_top = m.Slo_top_S(ij,1); % Slope at pixel ij
    Aspect = m.Aspect_S(ij,1); % Aspect at pixel ij
    clear SvF_S Ct_S Aspect_S

    cos_fst = cos(atan(Slo_top)).*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;

    clear zeta_S 
  
    SAB1(sin(h_S) <= 0.10)  =  0;
    SAB2(sin(h_S) <= 0.10)  =  0;
    SAD1(sin(h_S) <= 0.10)  =  0;
    SAD2(sin(h_S) <= 0.10)  =  0;
    PARB(sin(h_S) <= 0.10)  =  0;
    PARD(sin(h_S) <= 0.10)  =  0;
    SAD1 = SAD1.*SvF + Ct.*rho_g.*((SAB1./sin(h_S)).*cos_fst + (1-SvF).*SAD1);
    SAD2 = SAD2.*SvF + Ct.*rho_g.*((SAB2./sin(h_S)).*cos_fst + (1-SvF).*SAD2);
    PARD = PARD.*SvF + Ct.*rho_g.*((PARB./sin(h_S)).*cos_fst + (1-SvF).*PARD);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1 =(SAB1./sin(h_S)).*cos_fst.*ShF;
    SAB2 =(SAB2./sin(h_S)).*cos_fst.*ShF;
    PARB = (PARB./sin(h_S)).*cos_fst.*ShF;

    %correSITEions, temporary
    SAB1(SAB1<0)=0;SAB2(SAB2<0)=0;PARB(PARB<0)=0;PARD(PARD<0)=0;
    SAD1(SAD1<0)=0;SAD2(SAD2<0)=0;
    SAB1(isnan(SAB1)) = 0;SAB2(isnan(SAB2)) = 0;SAD1(isnan(SAD1)) = 0;SAD2(isnan(SAD2)) = 0;
    PARB(isnan(PARB)) = 0;PARD(isnan(PARD)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   SvF=1; %% Sky View factor = 1 if topography is not considered
   Slo_top_S = 0;
   Slo_top = 0;
end

SAB1(SAB1<0)=0;SAB2(SAB2<0)=0;PARB(PARB<0)=0;PARD(PARD<0)=0;
    SAD1(SAD1<0)=0;SAD2(SAD2<0)=0;
    SAB1(isnan(SAB1)) = 0;SAB2(isnan(SAB2)) = 0;SAD1(isnan(SAD1)) = 0;SAD2(isnan(SAD2)) = 0;
    PARB(isnan(PARB)) = 0;PARD(isnan(PARD)) = 0;

%% LAND COVER
deb_bool = glacier_points.deb_bool(point_id);
if deb_bool ==1
    ksv = 7;    %debris 
elseif deb_bool ==0
    ksv = 8;    %clean ice
end


%%% 1 Fir (evergr.)
%%% 2 Larch (decid.)
%%% 3 Grass C3
%%% 4 Shrub (decid.)
%%% 5 Broadleaf evergreen
%%% 6 Broadleaf deciduous
%%% 7 Rock
ksv(ksv==8)=7; %%% 8 Ice = 7 Rock

%land cover partition
cc_max = 1; %% one vegetation types

switch ksv
    case 1
        %%%% LAND COVER PARTITION  Fir - evergreen
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [1 0 0 0 0 0]>0;  
    case 2
        %%%% LAND COVER PARTITION  Larch - deciduous
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];  
        cc=length(Ccrown);%% Crown area
        II = [0 1 0 0 0 0]>0;  
    case 3
        %%%% LAND COVER PARTITION  Grass C3
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 1 0 0 0]>0;  
    case 4
        %%%% LAND COVER PARTITION  Shrub dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.1; Ccrown = [0.9];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 1 0 0]>0;  
    case 5
        %%%% LAND COVER PARTITION  broadleaf evergreen vegetation dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 0 1 0]>0;  
    case 6
        %%%% LAND COVER PARTITION  broadleaf deciduous vegetation dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 0 0 1]>0;  
    case 7
        %%%% LAND COVER PARTITION  Rock
        Cwat = 0; Curb = 0.0 ; Crock = 1.0;
        Cbare = 0.0; Ccrown = [0.0];
        cc=length(Ccrown);%% Crown area
        II = [ 0 0 0 1 0 0]>0;  


    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        return
end

zatm = max(zatm_surface(II)); %choose correct atmospheric reference height







%%% Initial snow depth

% Initial snow albedo
% if ~exist('SNOWALB','var')
%     SNOWALB = SNOWD;
%     SNOWALB(SNOWD>0) = 0.6;
% else 
%     SNOWALB=reshape(SNOWALB,num_cell,1);
% end 

%% INITIAL CONDITIONS AND PARAMETERS 
%%% SOIL
ms=10 ; %% 11 ; 
ms_max = 10; %% Number of soil layers
%%% DEBRIS
md_max = 10;
num_cell = 1; m_cell =1; n_cell=1;
MASKn = [1];
GLH = [100];
SNOWD = [1]; %%%% 1m is probably too thick!!!
SNOWALB = [0.6];
cellsize = 100;
Psan= 0.5;Pcla= 0.1;Porg = 0;
SOIL_TH = [200];
ij=1;
dbThick=glacier_points.deb_thickness_m(point_id)*1000;  %% [mm]

% if exist(out, 'file') == 2      %ADD BACK IN AFTER TESTING
% load(out);
% else
INIT_COND_v2(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Slo_top_S,ksv,Ca,SNOWD,SNOWALB,out);
load(out);
% end

%define parameter file
PARAM_IC = strcat(cd,'/Inputs/MOD_PARAM_point.m');





%% RUN MODEL
disp('starting_sim')
MAIN_FRAME; % Launch the main frame of T&C. Most of the things happen in this line of code
disp('ending_sim')






%% Output manager
Param_t = table(Lat,Lon,Zbas,dbThick,'VariableNames',{'Lat','Lon','Zbas','dbThick'});
Param_t = [Param_t, struct2table(SnowIce_Param), struct2table(Deb_Par)];
Param_t = rows2vars(Param_t);
Param_t = renamevars(Param_t,{'OriginalVariableNames','Var1'},{'Parameter','Value'});

%%post-compute sublimation from ESN
SSN = ESN.*(Ts<0);

% Here I manually choose the T&C outputs I want to save at each point. 

Outputs_t = table(Date,EICE,ESN,SND,SWE,...
Ta,Ws,U,N,SAD1+SAD2+SAB1+SAB2,Pre,Pr,Pr_sno,ALB,Smelt,Imelt,SSN,ICE,ET,ros,'VariableNames',{ ...
'Date','EICE','ESN','SND','SWE',...
'Ta','Ws','U','N','Rsw',...
'Pre','Pr','Pr_sno','Albedo','Smelt','Imelt','SSN','ICE','ET','ros'});

if output_daily == 1 % If I want to only save daily aggregated output 

Outputs_tt = table2timetable(Outputs_t);

Outputs_ds = retime(Outputs_tt,'daily',@nansum);
Outputs_dm = retime(Outputs_tt,'daily',@nanmean);

Outputs_d = Outputs_dm; 
Outputs_d.Pr = Outputs_ds.Pr;
Outputs_d.Pr_sno = Outputs_ds.Pr_sno;

Outputs_t = timetable2table(Outputs_d);

end 

parquetwrite( strcat(outlocation,'/',num2str(point_id),'_results.parquet'),Outputs_t)
writetable(Param_t, strcat(outlocation,'/',num2str(point_id),'_param.txt') )


%% Quick run evaluation

%%%%%% plot output 
% fi1 = figure('Renderer', 'painters', 'Position', [141.6667 244.3333 920.0000 370.6667]);

% plot(Outputs_t.Date,Outputs_t.SND)
% title(['ERA5L forcing + TA @ GP ' num2str(point_id)])
