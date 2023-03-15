%merge SeaFET raw_pH with TS data from SBE37:

cd 268_raw_data
load ob_268_raw_pH_
raw_pH = raw_pH_268; 
cd ..

%Time is UTC
year = floor(raw_pH(:,1)./1000);
jdy = raw_pH(:,1) - (year.*1000);
hr = floor(raw_pH(:,2));
mn = (raw_pH(:,2) - hr).*60;
pHtime = jul2date(jdy,year);
pHtime = datevec(pHtime);
pHtime = datenum(pHtime(:,1),pHtime(:,2),pHtime(:,3),hr,mn,ones(length(mn),1));

cd ..
cd CTD_raw_data
load TS_OB_2022_12_12.csv
ctd = TS_OB_2022_12_12; 
cd ..
cd 268_processing

%column headings:
%(1)date (excel SDN UTC)
%(2)Temp (deg C)
%(3)Oxygen (ml/l)
%(4)Salinity

date = excel2sdn(ctd(:,1)); % convert excel datenumber to matlab
ctd(:,1) = date;
xtick = [min(date):30:max(date)]; 

%% remove bad ctd salinity data

% remove all really low values
ck = find(ctd(:,4)< 27);
ctd(ck,4) = NaN;

% remove values from when cell was clogged and salinity was stable and low
% (Nov 17-Dec 3)
ck = find(ctd(:,1) > 738477.083 & ctd(:,1) < 738494); 
ctd(ck,4) = NaN;

% remove bad salinity data from March, cell probably clogged

ck = find(ctd(:,1) > 738594 & ctd(:,1) < 738606 & ctd(:,4)< 29); 
ctd(ck,4) = NaN;


%some spikes also removed manually, often times it would take salinity a
%while to get back to normal after sensor was pulled out of the water

%% View data

fig = figure;
set(fig,'DefaultAxesFontSize',15);

subplot(3,1,1)
title('Salinity')
plot(ctd(:,1),ctd(:,4),'b-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(date) max(date)]);
ylabel('Salinity')

subplot(3,1,2)
title('Temperature')
plot(ctd(:,1),ctd(:,2),'b-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(date) max(date)]);
ylabel('Temp (deg C)')

subplot(3,1,3)
title('Oxygen')
plot(ctd(:,1),ctd(:,3)./10,'b-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(date) max(date)]);
ylabel('Oxygen (mL/L)')

%% merge data with SeaFET

% (1)Datenumber	
% (2)Temp[C]	
% (3)Sal

% check temperatures match up
figure;
plot(pHtime,raw_pH(:,5),'k-')
hold on
plot(ctd(:,1),ctd(:,2),'r-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))

%interpolate CTD T and S to raw_pH timestep:

intT = interp1(ctd(:,1),ctd(:,2),pHtime);
intS = interp1(ctd(:,1),ctd(:,4),pHtime);
intOx = interp1(ctd(:,1),ctd(:,3),pHtime);

figure;
subplot(2,1,1)
plot(pHtime,raw_pH(:,5),'k-')
hold on
plot(pHtime,intT,'r-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,5))
subplot(2,1,2)
plot(pHtime,intT-raw_pH(:,5),'k.')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,5))

figure
plotregression(intT,raw_pH(:,5), 'Regression');
xlabel('CTD T (interpolated)')
ylabel('SeaFET T')

m = mean(intT-raw_pH(:,5),'omitnan')

%%
merged = [pHtime raw_pH(:,[ 2 3 4 6 7]) intT intS intOx];


%(1)Datetime (UTC, Matlab)
%(2)time decimal hour
%(3)7.80829, - pH int
%(4)7.78879, - pH ext
%(5)0.05975991, - isfet internal V
%(6)-0.85879952, - isfet external V
%(7)interpolated CTD Temp
%(8)interpolated CTD Salinity
%(9)interpolated CTD Oxygen

save merged_268_pH_CTD merged -ascii -double -tabs
