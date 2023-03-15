%Read raw SeaFET data
%Wiley Evans w/ auto-generated MATLAB code using import data button on
%toolbar

%Nov 12, 2015

cd 268_raw_data

%Full ASCII data file format:

%(1)SATPHA0001, - Header
%(2)2014083, - Date YYYYDDD
%(3)13.0932646, - Time decimal hour
%(4)7.80829, - pH int
%(5)7.78879, - pH ext
%(6)20.5221, - isfet thermistor
%(7)20.5038, - CTD T
%(8)32.6027, - CTD S
%(9)4.669, - CTD O
%(10)10.523, - CTD P
%(11)0.05975991, - isfet internal V
%(12)-0.85879952, - isfet external V
%(13)1.13056064, - isfet thermistor V
%(14)11.326, - supply V
%(15)114, - supply current
%(16)0.5, - electronics compartment relative humbidty
%(17)4.882, - internal 5V power supply
%(18)9.602, - main battery pack V
%(19)6.106, - internal isolated supply V
%(20)5.903, - isolated battery pack V
%(21)100, - substrate leakage current
%(22)100, - counter electrode leakage current
%(23)-0.95693927, - counter electrode V
%(24)0x0000, - status
%(25)7 - check sum

%-----------------------
%Elements to Update:
% - File name
filename = '*000.CSV*'; %adjust the word between the astericks to be something found in all file names

% - Project identifier
proj = 'ob_268_'; %used to save the combined data

% - Lines per file
lines = 1500; %populate this with a number greater than the maximum number of lines in a single file
%this is used to preallocate the space for the combined data set

% - Number of header lines in files
HL = 8;

delimiter = ',';
startRow = 2;

% - Number of columns (to be read from file)
num_col = 17; %25 for full fill, only extracting 17

% - Column format of files
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%Create list of files
files = dir(filename); 

%Process file by file
dataset = zeros((lines*length(files)), num_col); %Preallocate array 
rw = 1;
for M=1:length(files)
    disp(['beginning ' files(M).name]);
    fid = fopen(files(M).name);
    fname = files(M).name;
    dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    % Convert the contents of columns containing numeric strings to numbers.
    % Replace non-numeric strings with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    
    for col=[1,2,3,4,5,6,11,12,13,14,15,16,17,18,19,20,21,22,23,25]
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;
                
                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end
    % Split data into numeric and cell columns.
    rawNumericColumns = raw(:, [1,2,3,4,5,6,11,12,13,14,15,16,17,18,19,20,21,22,23,25]);
    rawCellColumns = raw(:, [7,8,9,10,24]);
    
    
    % Replace non-numeric cells with NaN
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
    
    % Allocate imported array to column variable names
    SATFHR = cell2mat(rawNumericColumns(:, 1));
    date = cell2mat(rawNumericColumns(:, 2));
    time = cell2mat(rawNumericColumns(:, 3)); %contains cal fac data from header
    VarName4 = cell2mat(rawNumericColumns(:, 4));
    VarName5 = cell2mat(rawNumericColumns(:, 5));
    VarName6 = cell2mat(rawNumericColumns(:, 6));
    VarName7 = rawCellColumns(:, 1);
    VarName8 = rawCellColumns(:, 2);
    VarName9 = rawCellColumns(:, 3);
    VarName10 = rawCellColumns(:, 4);
    VarName11 = cell2mat(rawNumericColumns(:, 7));
    VarName12 = cell2mat(rawNumericColumns(:, 8));
    VarName13 = cell2mat(rawNumericColumns(:, 9));
    VarName14 = cell2mat(rawNumericColumns(:, 10));
    VarName15 = cell2mat(rawNumericColumns(:, 11));
    VarName16 = cell2mat(rawNumericColumns(:, 12));
    VarName17 = cell2mat(rawNumericColumns(:, 13));
    VarName18 = cell2mat(rawNumericColumns(:, 14));
    VarName19 = cell2mat(rawNumericColumns(:, 15));
    VarName20 = cell2mat(rawNumericColumns(:, 16));
    VarName21 = cell2mat(rawNumericColumns(:, 17));
    VarName22 = cell2mat(rawNumericColumns(:, 18));
    VarName23 = cell2mat(rawNumericColumns(:, 19));
    VarName24 = rawCellColumns(:, 5);
    VarName25 = cell2mat(rawNumericColumns(:, 20));
    dataset(rw:(rw+length(date)-1),:) = [date, time, VarName4, VarName5, VarName6, VarName11, VarName12, VarName13, VarName14, VarName15, VarName16, VarName17, VarName18, VarName19, VarName20, VarName21, VarName22];
    rw = rw + length(date);
    fclose(fid);
end
dataset = dataset(1:(rw-1),:); %Remove excess rows
ck = isnan(dataset(:,1)); %Remove NaNs
ck = find(ck == 0);
dataset = dataset(ck,:);
% yr = yr(1:(row-1),1); %Remove excess rows
% dataset = [dataset,yr];
% dataset = sortrows(dataset,30); %intitial sort due to multiple years of data in dataset

raw_pH_268 = dataset; clear dataset;

%test plot of data:
%plot1 is internal pH
%plot2 is supply V
%plot3 is leak electronics humidity

year = floor(raw_pH_268(:,1)./1000);
jdy = raw_pH_268(:,1) - (year.*1000);
hr = floor(raw_pH_268(:,2));
mn = (raw_pH_268(:,2) - hr).*60;

pHtime = jul2date(jdy,year);
pHtime = datevec(pHtime);
pHtime = datenum(pHtime(:,1),pHtime(:,2),pHtime(:,3),hr,mn,ones(length(mn),1));

figure; 
plot(pHtime,raw_pH_268(:,3),'k-')
xlabel('datenum')
ylabel('pH int')

figure; 
plot(pHtime,raw_pH_268(:,9),'k-')
xlabel('datenum')
ylabel('supply V')

figure; 
plot(pHtime,raw_pH_268(:,11),'k-')
xlabel('datenum')
ylabel('electronics compartment relative humidity')

%Save data as matlab file
saveas = strcat(proj, 'raw_pH_');
save(saveas,'raw_pH_268', '-v7.3');

%Clear extraneous variables (keeping these variables in workspace may be
%useful as a troubleshooting tool)
clear AA; clear dataset; clear fid; clear filename; clear files; clear format; clear HL; clear M; clear row;
clear lines; clear ans; clear num_col; clear raw_xCO2; clear saveas;
clear proj; clear year; clear VarName4; clear VarName5; clear VarName6; clear VarName7; clear VarName8;
clear VarName9; clear VarName10; clear VarName11; clear VarName12; clear VarName13; clear VarName14; clear VarName15;       
clear VarName16; clear VarName17; clear VarName18; clear VarName19; clear VarName20; clear VarName21;      
clear VarName22; clear VarName23; clear VarName24; clear VarName25;  
clear ck; clear col; clear dataArray; clear date; clear delimiter; clear fname; clear formatSpec; clear hr
clear invalidThousandsSeparator; clear jdy; clear me; clear mn; clear numbers; clear numericData; clear pHtime
clear R; clear raw; clear raw_pH; clear rawCellColumns; clear rawData; clear rawNumericColumns; clear regexstr;
clear result; clear rw; clear SATFHR; clear startRow; clear time

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