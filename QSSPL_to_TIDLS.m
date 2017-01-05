%{
MIT License

Copyright (c) [2016] [Mallory Ann Jensen, jensenma@alum.mit.edu]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}

%This function reads the QSSPL data (raw) for a certain sample, where lifetime
%measurements have been taken and processed separately.
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160811\A23-6';
%Filenames in order of temperature
filenames = {[dirname '\A23-6_25C.txt'];...
    [dirname '\A23-6_100C.txt'];...
    [dirname '\A23-6_150C.txt'];...
    [dirname '\A23-6_200C.txt']};

%!!! This is a placeholder
sample_param_filename = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160804\For processing\python_inputs.mat';
info = load(sample_param_filename); 

%For each temperature, load the data
for i = 1:length(filenames)
    format_for_TIDLS(filenames{i},sample_param_filename,dirname);
    [tau_mean,deltanq]=average_tau([dirname '\Raw_data.mat'],dirname);
end

%% Now try averaging for a series 
clear all; close all; 
dirname = 'C:\Users\Malloryj\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160905\16-6-28-P-2';
sample = '16-6-28-P-2'; 
Ctonum = '_IRap_1-1_';
time_before = .001; %s
T = [25]; 
num_repeats = [5]; 
[dataSave] = average_QSSPL(dirname,sample,T,num_repeats,time_before,Ctonum);

%Now, given this data, we want to write this to a new file which can be
%read into QSSPL
for i = 1:length(T)
    matrix_to_write = dataSave{i}; 
    write_QSSPL_analyser(dirname,sample,T(i),matrix_to_write,Ctonum);
end
    
%% Now we prepare the data in our expected TIDLS format
%Tell us how the files were entered in terms of temperature. These are the
%actual sample temperatures in K. 
direction = 'ascending'; 
sample_name = '22-102-8'; 
T = [466 418 370 298]; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\22-102-8\flash only';

%Read and save the raw data
data_filename = [dirname '\' sample_name '_all_' direction '.txt']; 
format_for_TIDLS(data_filename,dirname,'PL');

%Now that we've saved the raw data, we need to save the sample parameters
%as expected. For this we read the .inf files
% if direction == 'ascending'
%     T_list = flip(T); 
% else
%     T_list = T; 
% end
% for i = 1:length(T)
%     %Make the .inf filename
%     inf_name = [dirname '\' sample_name '_' num2str(T) 'C_averaged.inf']; 
%     %Read the inf

%For the moment, just manually enter the information
for i = 1:length(T)
    temperature{i} = T(i)-273.15; 
    resistivity{i} = 1000; 
    thickness{i} = .0173; 
    OC{i} = 0.8; 
    calibration_factor{i} = 1; 
    doping{i} = 8.5e15; 
    fileListShort{i} = [sample_name '_' num2str(T(i)) 'C'];
end
info = struct('filename',fileListShort,'thickness',thickness,'resistivity',...
    resistivity,'measured_resistivity',resistivity,'optical_constant',OC,...
    'calibration_factor',calibration_factor,'temperature',temperature,'doping',doping);
save([dirname '\meas_info.mat'],'info'); 

%% Stitch together data from multiple measurements
clear all; close all; 
filename = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\22-25-8\25C\Attempt 2\22-25-8.txt'; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\22-25-8\25C\Attempt 2';
% format_for_TIDLS(filename,dirname,'PC');
load([dirname '\PC\Raw_data.mat']); 
dataPC = dataSave; 
load([dirname '\PL\Raw_data.mat']);
dataPL = dataSave; 
figure;
for i = 1:length(dataPC)
    to_plot = dataPC{i}; 
    loglog(to_plot(:,1),to_plot(:,2),'-','LineWidth',2); 
    hold all; 
    to_plot = dataPL{i}; 
    loglog(to_plot(:,1),to_plot(:,2),'--','LineWidth',2); 
    hold all; 
end
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20); 
%Now try to combine all the data. Just try a simple combination
for i = 1:length(dataPC)
    if i == 1
        togetherPC = dataPC{i}; 
        togetherPL = dataPL{i}; 
    else
        togetherPC = [togetherPC;dataPC{i}];
        togetherPL = [togetherPL;dataPL{i}];
    end
end
togetherPC = sortrows(togetherPC,1); 
togetherPL = sortrows(togetherPL,1); 
hold all; 
loglog(togetherPC(:,1),togetherPC(:,2),'k-','LineWidth',2); 
hold all; 
loglog(togetherPL(:,1),togetherPL(:,2),'k--','LineWidth',2); 

%% 
clear all; close all; 
sample = 'A41-6'; 
dirname = 'C:\Users\Malloryj\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\A41-6\Summary files'; 
T = [-100 -50 25 100]; %C
dataStore = cell(length(T),1); 
figure; 
for i = 1:length(T)
%     filename = [dirname '\' num2str(T(i)) 'C\' sample '.txt']; 
%     filename = [dirname '\' num2str(T(i)) 'C\' sample '_' num2str(T(i)) 'C.txt'];
    filename = [dirname '\' num2str(T(i)) 'C\' sample '_' num2str(T(i)) 'C_PC.txt'];
    format_for_TIDLS(filename,[dirname '\' num2str(T(i)) 'C'],'PC');
    load([dirname '\' num2str(T(i)) 'C\Raw_data.mat']); 
    data_now = dataSave; 
    for j = 1:length(data_now)
        if j == 1
            togetherPC = data_now{j}; 
        else
            togetherPC = [togetherPC;data_now{j}];
        end
    end
    togetherPC = sortrows(togetherPC,1); 
    loglog(togetherPC(:,1),togetherPC(:,2),'k-','LineWidth',2); 
    hold all; 
    dataStore{i} = togetherPC; 
end

T = [177 225 296 368]; 
for i = 1:length(T)
    temperature{i} = T(i)-273.15; 
    resistivity{i} = 0.95; 
    thickness{i} = .0193; 
    OC{i} = 0.7; 
    calibration_factor{i} = 1; 
    doping{i} = 1.6e16; 
    fileListShort{i} = [sample '_' num2str(T(i)) 'C'];
end
info = struct('filename',fileListShort,'thickness',thickness,'resistivity',...
    resistivity,'measured_resistivity',resistivity,'optical_constant',OC,...
    'calibration_factor',calibration_factor,'temperature',temperature,'doping',doping);
save([dirname '\meas_info.mat'],'info'); 

dataSave = dataStore; 
save([dirname '\Raw_data.mat'],'dataSave'); 

