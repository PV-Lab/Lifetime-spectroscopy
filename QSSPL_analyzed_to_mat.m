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

%% Format data analysed by QSSPL-analyser into something that can be analyzed by TIDLS
%Ideally the data for all relevant temperatures has already been analyzed
%and can now be constructed together into one .mat file

clear all; close all; clc; 

%Sample name
sample = '22-102-8'; 
%Where the files are stored, after this directory there should be folders
%for each temperature
dirname = 'C:\Users\Mallory Jensen\Dropbox (MIT)\TIDLS data - UROP\Data\sorted by sample\22 102 8\analyzed data\PL'; 
%Temperatures we're looking at
T = [-75 -25 25 100]; %C
%Which data do we want, PC or PL
PCorPL = 'PC'; 

%Now the actual sample values
%Actual sample temperatures, used for analysis and in units of Kelvin
sampleT = [298.15 367]; 
thick = .028; %thickness, cm
res = 1.1; %resistivity, ohm-cm
opt = 0.9; %This is the 1 minus the reflectivity value used for flash measurements
N_D = 1.5e15; %doping level, n-type, cm^-3

dataStore = cell(length(T),1); 
figure; 
for i = 1:length(T)
    try
        filename = [dirname '\' num2str(T(i)) 'C\' sample '_' num2str(T(i)) 'C.txt']; 
        format_for_TIDLS(filename,[dirname '\' num2str(T(i)) 'C'],PCorPL);
        load([dirname '\' num2str(T(i)) 'C\Raw_data.mat']); 
        data_now = dataSave; 
        for j = 1:length(data_now)
            if j == 1
                together = data_now{j}; 
            else
                together = [together;data_now{j}];
            end
        end
        together = sortrows(together,1); 
        loglog(together(:,1),together(:,2),'k-','LineWidth',2); 
        hold all; 
        dataStore{i} = together; 
    catch
        disp(['Could not load data for T = ' num2str(T(i)) 'C']);
    end
end

for i = 1:length(sampleT)
    temperature{i} = sampleT(i)-273.15; 
    resistivity{i} = res; 
    thickness{i} = thick; 
    OC{i} = opt; 
    calibration_factor{i} = 1; %this doesn't matter here but we need it to prevent errors
    doping{i} = N_D; 
    fileListShort{i} = [sample '_' num2str(sampleT(i)) 'C'];
end
info = struct('filename',fileListShort,'thickness',thickness,'resistivity',...
    resistivity,'measured_resistivity',resistivity,'optical_constant',OC,...
    'calibration_factor',calibration_factor,'temperature',temperature,'doping',doping);
save([dirname '\meas_info.mat'],'info'); 

dataSave = dataStore; 
save([dirname '\Raw_data.mat'],'dataSave'); 
