% plot diagnostics & calculate in situ coefficients

cd 268_raw_data
load ob_268_raw_pH_
raw_pH = raw_pH_268; 
cd ..
%column headings:
%(1)date YYYYDDD
%(2)time decimal hour
%(3)7.80829, - pH int
%(4)7.78879, - pH ext
%(5)20.5221, - isfet thermistor
%(6)0.05975991, - isfet internal V
%(7)-0.85879952, - isfet external V
%(8)1.13056064, - isfet thermistor V
%(9)11.326, - supply V
%(10)114, - supply current
%(11)0.5, - electronics compartment relative humbidty
%(12)4.882, - internal 5V power supply
%(13)9.602, - main battery pack V
%(14)6.106, - internal isolated supply V
%(15)5.903, - isolated battery pack V
%(16)100, - substrate leakage current
%(17)100, - counter electrode leakage current

%Time is UTC
year = floor(raw_pH(:,1)./1000);
jdy = raw_pH(:,1) - (year.*1000);
hr = floor(raw_pH(:,2));
mn = (raw_pH(:,2) - hr).*60;
pHtime = jul2date(jdy,year);
pHtime = datevec(pHtime);
pHtime = datenum(pHtime(:,1),pHtime(:,2),pHtime(:,3),hr,mn,ones(length(mn),1));

% 1.Diagnostic plots

figure; 
plot(pHtime,raw_pH(:,3),'k-')
xlabel('datenum')
hold on

figure; 
plot(pHtime,raw_pH(:,9),'k-')
xlabel('datenum')
ylabel('supply V')
xtick = [min(pHtime):10:max(pHtime)];
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);

figure; 
plot(pHtime,raw_pH(:,11),'k-')
xlabel('datenum')
ylabel('electronics compartment relative humidity')
xtick = [min(pHtime):10:max(pHtime)];
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
set(gca,'xlim',[min(pHtime) max(pHtime)]);

figure;
plot(pHtime,raw_pH(:,6),'b-')
hold on
plot(pHtime,raw_pH(:,7),'g-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
ylabel('Voltage')
legend('pH_i_n_t V','pH_e_x_t V');
set(gca,'xlim',[min(pHtime) max(pHtime)]);

figure;
plot(pHtime,raw_pH(:,3),'b-')
hold on
plot(pHtime,raw_pH(:,4),'g-')
set(gca,'xtick',xtick,'xticklabel',datestr(xtick,6))
ylabel('pH')
legend('pH_i_n_t','pH_e_x_t');
set(gca,'xlim',[min(pHtime) max(pHtime)]);

%% Calculate in situ calibration coefficient: Segment 1
% 
% in situ calibration using second sampling event June 10th 2021

%%%%%Adjust according to deployment:
%cal data point time:
calD = datenum(2021,6,10,23,30,00);   % UTC
calTemp = 11.5;  %probe temp
calpH = 7.786626037;  %discrete sample pHT from calibration data file (used middle triplicate)
calS = 29.0;   % YSI S
k2_i = [-1.101000e-03]; %get from header of raw data file (CAL_PHINT_SLOPE_COEFF)
k2_e = [-1.048000e-03]; %get from header of raw data file (CAL_PHEXT_SLOPE_COEFF)

%first - internal V
figure; 
plot([calD calD],[min(raw_pH(:,6)) max(raw_pH(:,6))],'r-','linewidth',4)
hold on
plot(pHtime,raw_pH(:,6),'ko','markerfacecolor','k','markersize',5)
xlabel('datenum')
ylabel('internal V')
%axis([calD-0.5 calD+10 -inf inf])
dtick = [calD-0.5:2/24:calD+0.5];
set(gca,'xtick',dtick,'xticklabel',datestr(dtick,0))
%pbaspect([7 1 1])
grid on
box on

%isolate voltage nearest cal sample:
lowCalD = calD-(5/1440);    %within +/- 5 minutes
highCalD = calD+(5/1440);

ck = find(pHtime > lowCalD & pHtime < highCalD);
VforCal = raw_pH(ck,:);
VforCalT = pHtime(ck);

%note: the number of rows here is the number of frames in a burst. Current
%settings for our sensor are average of 10 measurements = 1 frame for burst
%of 30 frames @ a rate of 1 frame/sec

%average burst of frames and use these voltages to calculate cal offsets
%(ko_i, ko_e)

VforCal = mean(VforCal,1);  %1x17
VforCalT = mean(VforCalT,1);

%Equations below are from SeaFET Manual 2.0 pg.12-13
%See Cale et al excel calculator for comparison
%calculate ko_i:
Snerst = (8.314472*(calTemp+273.15)*log(10))/96485.3415;
ko_i = VforCal(6)-((calpH*Snerst)+(k2_i*(calTemp+273.15)))

%calculate ko_e:
Cl = (0.99889/35.453)*(calS/1.80655); %total chloride
Adh = (0.00000343*(calTemp^2)) + (0.00067524*calTemp) + 0.49172143; %Debye-Huckel constant for activity of HCl
I = (19.924*calS)/(1000-1.005*calS); %Ionic strength
St = (0.1400/96.062)*(calS/1.80655); %total sulfate
Ks = (1-0.001005*calS)*exp((-4276.1/(calTemp+273.15))+141.328-(23.093*log((calTemp+273.15)))...
    +(-13856/(calTemp+273.15) + 324.57 - 47.986*log((calTemp+273.15)))*sqrt(I)...
    +((35474/(calTemp+273.15)-771.54+114.723*log(calTemp+273.15))*I)-((2698/(calTemp+273.15))*(I^1.5))...
    +((1776/(calTemp+273.15))*(I^2))); %acid dissociation constant of HSO4-
HCL = ((-Adh*sqrt(I))/(1+1.394*sqrt(I)))+(0.08885-0.000111*calTemp)*I; %HCL = log10(HCL); %logarithm of HCL activity coefficient

ko_e = VforCal(7)-((calpH+(log10(1+(St/Ks))-(2*HCL)-log10(Cl)))*Snerst + k2_e*(calTemp+273.15))

cal_data_1 = [ko_i k2_i ko_e k2_e];

%% Calculate in situ calibration coefficient: Segment 2
% 
% in situ calibration using second sampling event: May 6th 2022

%%%%%Adjust according to deployment:
%cal data point time:
calD = datenum(2022,5,6,18,30,00);   % UTC
calTemp = 10.2;  %probe temp
calpH = 7.829859283;  %discrete sample pHT from calibration data file (used middle triplicate)
calS = 29.4;   % YSI S
k2_i = [-1.101000e-03]; %get from header of raw data file (CAL_PHINT_SLOPE_COEFF)
k2_e = [-1.048000e-03]; %get from header of raw data file (CAL_PHEXT_SLOPE_COEFF)

%first - internal V
figure; 
plot([calD calD],[min(raw_pH(:,6)) max(raw_pH(:,6))],'r-','linewidth',4)
hold on
plot(pHtime,raw_pH(:,6),'ko','markerfacecolor','k','markersize',5)
xlabel('datenum')
ylabel('internal V')
%axis([calD-0.5 calD+10 -inf inf])
dtick = [calD-0.5:2/24:calD+0.5];
set(gca,'xtick',dtick,'xticklabel',datestr(dtick,0))
%pbaspect([7 1 1])
grid on
box on

%isolate voltage nearest cal sample:
lowCalD = calD-(5/1440);    %within +/- 5 minutes
highCalD = calD+(5/1440);

ck = find(pHtime > lowCalD & pHtime < highCalD);
VforCal = raw_pH(ck,:);
VforCalT = pHtime(ck);

%note: the number of rows here is the number of frames in a burst. Current
%settings for our sensor are average of 10 measurements = 1 frame for burst
%of 30 frames @ a rate of 1 frame/sec

%average burst of frames and use these voltages to calculate cal offsets
%(ko_i, ko_e)

VforCal = mean(VforCal,1);  %1x17
VforCalT = mean(VforCalT,1);

%Equations below are from SeaFET Manual 2.0 pg.12-13
%See Cale et al excel calculator for comparison
%calculate ko_i:
Snerst = (8.314472*(calTemp+273.15)*log(10))/96485.3415;
ko_i = VforCal(6)-((calpH*Snerst)+(k2_i*(calTemp+273.15)))

%calculate ko_e:
Cl = (0.99889/35.453)*(calS/1.80655); %total chloride
Adh = (0.00000343*(calTemp^2)) + (0.00067524*calTemp) + 0.49172143; %Debye-Huckel constant for activity of HCl
I = (19.924*calS)/(1000-1.005*calS); %Ionic strength
St = (0.1400/96.062)*(calS/1.80655); %total sulfate
Ks = (1-0.001005*calS)*exp((-4276.1/(calTemp+273.15))+141.328-(23.093*log((calTemp+273.15)))...
    +(-13856/(calTemp+273.15) + 324.57 - 47.986*log((calTemp+273.15)))*sqrt(I)...
    +((35474/(calTemp+273.15)-771.54+114.723*log(calTemp+273.15))*I)-((2698/(calTemp+273.15))*(I^1.5))...
    +((1776/(calTemp+273.15))*(I^2))); %acid dissociation constant of HSO4-
HCL = ((-Adh*sqrt(I))/(1+1.394*sqrt(I)))+(0.08885-0.000111*calTemp)*I; %HCL = log10(HCL); %logarithm of HCL activity coefficient

ko_e = VforCal(7)-((calpH+(log10(1+(St/Ks))-(2*HCL)-log10(Cl)))*Snerst + k2_e*(calTemp+273.15))

cal_data_2 = [ko_i k2_i ko_e k2_e];


