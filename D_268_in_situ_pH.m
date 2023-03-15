%re-compute pH using CTD data and in situ calibration coefficients

data = load('merged_268_pH_CTD');

%(1)Datetime (UTC, Matlab)
%(2)time decimal hour
%(3)7.80829, - pH int
%(4)7.78879, - pH ext
%(5)0.05975991, - isfet internal V
%(6)-0.85879952, - isfet external V
%(7)interpolated CTD Temp
%(8)interpolated CTD Salinity
%(9)interpolated CTD Oxygen

pHtime = data(:,1);

%% %average 30 burst samples:
dum = diff(pHtime);
figure
plot(dum)
ck = find(dum > 0.003); % 5 minute sampling interval (5 mins = 0.0035)
n = length(pHtime);
ck = [0;ck;n];
%prealocate the array:
mdata = zeros(length(ck)-1, size(data,2));
for i = 2:length(ck)
    insert = nanmean(data(ck(i-1)+1:ck(i),:));
    mdata(i-1,:) = insert;
end

%replace data file and recompute time:
data = mdata;
pHtime = data(:,1);

%% clean up SeaFET data

figure;
xtick = [min(data(:,1)):15:max(data(:,1))]; 
plot(data(:,1),data(:,3),'b-')
hold on
plot(data(:,1),data(:,4),'g-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
ylabel('pH')
legend('pH_i_n_t','pH_e_x_t');
set(gca,'xlim',[min(pHtime) max(pHtime)]);

% remove bad data associated with service visits (pH and temp spikes, dates confirmed to be service visits)

ck = find(data(:,1) > 738377.83 & data(:,1) < 738378.00);
data(ck,[3 4 5 6]) = NaN;

ck1 = find(data(:,1) > 738415.73 & data(:,1) < 738415.83);
data(ck1,[3 4 5 6]) = NaN;

ck2 = find(data(:,1) > 738486.75 & data(:,1) < 738486.85);
data(ck2,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738493.73 & data(:,1) < 738493.83);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738486.7084 & data(:,1) < 738486.70856);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738493.72939 & data(:,1) < 738493.8335);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738575.750227232 & data(:,1) < 738575.771060458);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738575.750227232 & data(:,1) < 738575.771060458);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738577 & data(:,1) < 738610);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738426.30 & data(:,1) < 738426.7293);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738486.854393931 & data(:,1) < 738487.14606054);
data(ck3,[3 4 5 6]) = NaN;

ck3 = find(data(:,1) > 738649.7502271396 & data(:,1) < 738649.79189);
data(ck3,[3 4 5 6]) = NaN;

%%
figure;
subplot(2,1,1)
plot(pHtime,data(:,3),'b-')
xtick = [min(pHtime):30:max(pHtime)];
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);
ylabel('V int')
pbaspect([4 1 1])
subplot(2,1,2)
plot(pHtime,data(:,4),'r-')
xtick = [min(pHtime):30:max(pHtime)];
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);
ylabel('V ext')
pbaspect([4 1 1])

figure;
subplot(2,1,1)
plot(pHtime,data(:,7),'b-')
xtick = [min(pHtime):30:max(pHtime)];
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);
hold on
plot(pHtime,data(:,9),'k-')
ylabel('T')
legend('SBE T [interpolated]','SBE O2 [interpolated]')
pbaspect([4 1 1])

subplot(2,1,2)
plot(pHtime,data(:,8),'r-')
xtick = [min(pHtime):30:max(pHtime)];
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);
ylabel('S')
pbaspect([4 1 1])

%% First, plot factory calibrated pH with bottle data to determined cal coefficient to use

cd ..
btl = load('OB_btl_data.csv');
cd 268_processing

%columns
%1 Date & time Collected (excel number) UTC
%2 In situ Temp
%3 YSI Sal
%4 CRM corrected TCO2
%5 correct TA
%6 pCO2@SST
%7 pH_T
%8 station (5= sensor, 1=kelp, 6=outside in channel, 3= owen bay, 7=creek)

btl_date = excel2sdn(btl(:,1));
btl(:,1) = btl_date;

ck = find(btl(:,8) == 5);
sensor = btl(ck,:);

ck = find(btl(:,8) == 1);
kelp = btl(ck,:);

ck = find(btl(:,8) == 6);
channel = btl(ck,:);

ck = find(btl(:,8) == 3);
ob = btl(ck,:);

ck = find(btl(:,8) == 7);
fw = btl(ck,:);

% plot electrode pH (factory calibrated) with bottle pH
figure;
plot(pHtime,data(:,4),'k-')
hold on
plot(pHtime,data(:,3),'b-')
hold on
plot(sensor(:,1),sensor(:,7),'ko','markerfacecolor','m','markersize',5);
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);
ylabel('pH')
legend('pH_e_x_t fact','pH_i_n_t fact','bottle pH');


% lots of bottle samples not shown as they don't have concurrent data from
% sensor because it was removed in QC due to spikes in pH or temp after
% being out of the water :(

% First segment: May 2021 - Feb 2022
% only good samples for calibration are the ones collected by Kyle on 
% June 10th 2021, visual comparison with rest show that internal electrode
% was pretty close to bottle pH for the whole segment

% Second segment: March 2022 - June 2022
% batteries changed followed by one month of data loss, calibration
% coefficients show reversal in correction for the internal electrode so
% need to do a seperate calibration

%% compute differences between discrete pH and SeaFET pH:
sensor_pH = interp1(data(:,1),data(:,3),sensor(:,1));  %%%%%%%%%%%%%%%%%%%%%%%3 = external, 4 = internal
diff_from_sensor_pH = sensor(:,7) - sensor_pH;
figure;
subplot(2,1,1)
plot(sensor(:,1),diff_from_sensor_pH,'bo','markerfacecolor','b','markersize',5);
hold on
x = [min(sensor(:,1)),max(sensor(:,1))];
y = [0,0];
plot(x,y,'k--','linewidth',2);
hold on
ylabel('btl pH - sensor pH')
set(gca,'ylim',[-0.12 0.12]);
set(gca,'xlim',[min(pHtime) max(pHtime)]);
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,3))
title('internal');

sensor_pH = interp1(data(:,1),data(:,4),sensor(:,1));  %%%%%%%%%%%%%%%%%%%%%%%3 = external, 4 = internal
diff_from_sensor_pH = sensor(:,7) - sensor_pH;
subplot(2,1,2)
plot(sensor(:,1),diff_from_sensor_pH,'bo','markerfacecolor','b','markersize',5);
hold on
x = [min(sensor(:,1)),max(sensor(:,1))];
y = [0,0];
plot(x,y,'k--','linewidth',2);
hold on
ylabel('btl pH - sensor pH')
set(gca,'ylim',[-0.12 0.12]);
set(gca,'xlim',[min(pHtime) max(pHtime)]);
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,3))
title('external');


%% Split data file into two segments

% segment 1
start = min(data(:,1));
fin = 738581; % March 1st 2022
ck1 = find(data(:,1) > start & data(:,1) < fin);
data1 = data(ck1,:);

% segment 2
ck2 = find(data(:,1) > fin);
data2 = data(ck2,:);

%% Calculate coefficients

% go back to step B to calculate cal coefficients as data has already been
% averaged at this point

cal_data_1 =[-1.08867467309650,-0.00110100000000000,-1.09306485566048,-0.00104800000000000]; % June 10th triplicate
cal_data_2 =[-1.09021810603268,-0.00110100000000000,-1.09258560418988,-0.00104800000000000]; % May 6th triplicate

%% Segment 1: calculate pH_int

% use in situ caldata and do TS correction
cal_data = cal_data_1;
Snerst = ((data1(:,7)+273.15).*(8.314472*log(10)))/96485.3415;
pH_int1 = (data1(:,5)- cal_data(1) - ((data1(:,7)+273.15).*cal_data(2)))./Snerst;

figure;
plot(data1(:,1),data1(:,3),'k-')
hold on
plot(data1(:,1),pH_int1,'r-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,5))
ylabel('pH_T')
legend('factory cal','in situ cal')
title('internal electrode')

%% Segment 1: calculate pH_ext

%calculate pH_ext:
Cl = (data1(:,8)./1.80655).*(0.99889/35.453); %total chloride
Adh = ((data1(:,7).^2).*0.00000343) + (data1(:,7).*0.00067524) + 0.49172143; %Debye-Huckel constant for activity of HCl
I = (data1(:,8).*19.924)./(1000-(data1(:,8).*1.005)); %Ionic strength
St = (data1(:,8)./1.80655).*(0.1400/96.062); %total sulfate
Ks = (1-(data1(:,8).*0.001005)).*exp((-4276.1./(data1(:,7)+273.15))+141.328-(23.093.*log((data1(:,7)+273.15)))...
    +(-13856./(data1(:,7)+273.15) + 324.57 - 47.986.*log((data1(:,7)+273.15))).*sqrt(I)...
    +((35474./(data1(:,7)+273.15)-771.54+114.723.*log(data1(:,7)+273.15)).*I)-((2698./(data1(:,7)+273.15)).*(I.^1.5))...
    +((1776./(data1(:,7)+273.15)).*(I.^2))); %acid dissociation constant of HSO4-
HCL = ((-Adh.*sqrt(I))./(1+1.394.*sqrt(I)))+(0.08885-0.000111.*data1(:,7)).*I; %HCL = log10(HCL); %logarithm of HCL activity coefficient

pH_ext1 = ((data1(:,6)-cal_data(3)-(cal_data(4).*(data1(:,7)+273.15)))./Snerst)+log10(Cl)+(2.*HCL)-(log10(1+(St./Ks)));

data1 = [data1 pH_int1 pH_ext1];

figure;
plot(data1(:,1),data1(:,4),'k-')
hold on
plot(data1(:,1),pH_ext1,'r-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,5))
ylabel('pH_T')
legend('factory cal','in situ cal')
title('external electrode')

%% Segment 2: calculate pH_int

% use in situ caldata and do TS correction
cal_data = cal_data_2;
Snerst = ((data2(:,7)+273.15).*(8.314472*log(10)))/96485.3415;
pH_int2 = (data2(:,5)- cal_data(1) - ((data2(:,7)+273.15).*cal_data(2)))./Snerst;

figure;
plot(data2(:,1),data2(:,3),'k-')
hold on
plot(data2(:,1),pH_int2,'r-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,5))
ylabel('pH_T')
legend('factory cal','in situ cal')
title('internal electrode')

%% Segment 2: calculate pH_ext

%calculate pH_ext:
Cl = (data2(:,8)./1.80655).*(0.99889/35.453); %total chloride
Adh = ((data2(:,7).^2).*0.00000343) + (data2(:,7).*0.00067524) + 0.49172143; %Debye-Huckel constant for activity of HCl
I = (data2(:,8).*19.924)./(1000-(data2(:,8).*1.005)); %Ionic strength
St = (data2(:,8)./1.80655).*(0.1400/96.062); %total sulfate
Ks = (1-(data2(:,8).*0.001005)).*exp((-4276.1./(data2(:,7)+273.15))+141.328-(23.093.*log((data2(:,7)+273.15)))...
    +(-13856./(data2(:,7)+273.15) + 324.57 - 47.986.*log((data2(:,7)+273.15))).*sqrt(I)...
    +((35474./(data2(:,7)+273.15)-771.54+114.723.*log(data2(:,7)+273.15)).*I)-((2698./(data2(:,7)+273.15)).*(I.^1.5))...
    +((1776./(data2(:,7)+273.15)).*(I.^2))); %acid dissociation constant of HSO4-
HCL = ((-Adh.*sqrt(I))./(1+1.394.*sqrt(I)))+(0.08885-0.000111.*data2(:,7)).*I; %HCL = log10(HCL); %logarithm of HCL activity coefficient

pH_ext2 = ((data2(:,6)-cal_data(3)-(cal_data(4).*(data2(:,7)+273.15)))./Snerst)+log10(Cl)+(2.*HCL)-(log10(1+(St./Ks)));

data2 = [data2 pH_int2 pH_ext2];

figure;
plot(data2(:,1),data2(:,4),'k-')
hold on
plot(data2(:,1),pH_ext2,'r-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,5))
ylabel('pH_T')
legend('factory cal','in situ cal')
title('external electrode')

%% Re-combine segments

data = vertcat(data1, data2);

%% Plot processed data with bottle samples & factory calibrated data

% plot 268 data raw vs. processed data
figure;
plot(data(:,1),data(:,4),'k-')
hold on
plot(data(:,1),data(:,11),'r-')
hold on
plot(data(:,1),data(:,3),'b-')
hold on
plot(data(:,1),data(:,10),'g-')
hold on
plot(sensor(:,1),sensor(:,7),'ko','markerfacecolor','m','markersize',5);
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
% set(gca,'ylim',[7.9 8.4]);
ylabel('pH')
legend('pH_e_x_t raw','pH_e_x_t processed','pH_i_n_t raw','pH_i_n_t processed','bottle pH');


%% Save data as matlab file
fet_data = data;
save('proc_pH_268','fet_data', '-v7.3');

%(1)datenumber (utc)
%(2)time decimal hour (utc)
%(3)7.80829, - pH int (factory cal)
%(4)7.78879, - pH ext (factory cal)
%(5)0.05975991, - isfet internal V
%(6)-0.85879952, - isfet external V
%(7)interpolated CTD T
%(8)interpolated CTD S
%(9)interpolated CTD O2
%(10)pH int (in situ cal 2021/22, TS corrected)
%(11)pH ext (in situ cal 2021/22, TS corrected)

%% Save data as csv file for Ondine
% exceldate = sdn2excel(fet_data(:,1));
% fet_data(:,1) = exceldate;
% 
% dlmwrite('proc_pH_454.csv',fet_data,'precision','%.10f') 
% save proc_pH_454 fet_data -ascii -double -tabs

