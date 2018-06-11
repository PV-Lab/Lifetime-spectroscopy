%{

MIT License

Copyright (c) [2018] [Mallory Ann Jensen, jensenma@alum.mit.edu & Ryan Sander, rmsander@mit.edu]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this 0permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}
clear all; close all; clc;

%This function is used to gather cell degradation data from Suns Voc
%Lifetime measurements.
%___________
%USER INPUTS
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Quasi Fermi level\UNSW collaboration\Cells at MIT\Suns-Voc\PERC 3-8';
base_sample_name = [directory '\PERC3-8_']; %E.g. 'BSF3-8'; don't include seconds
degradation_seconds = [0,500,1000,2500,5000,7500,7501,10000,23510,23511,39170,68090,68091,98330,98331,130520,160310,245570,269030,332570,352910,416390,446960,524660]; %Enter as a list where each index is a time measurement (e.g. 500 corresponds to a lifetime measurement at 500 seconds)
% degradation_seconds = [0 5 10 20 40 70 110 160 220 300 400 530 700 701 900 1150 1400 1650 1900 2150 2400 2700 3000 3250 3500 3800 4100 4400];
% degradation_seconds = [0 5 10 20 40 70 110 160 220 300 400 530 700 701 900 1400 1650 1900 2150 2400 2700 3000 3250 3500 3800 4100 4400];
deltan_target = 6e14; 
%___________ABOVE ARE THE ONLY TWO USER INPUTS
%Initialize output vectors and tensors
num_files = length(degradation_seconds);
VOCs = zeros(num_files,1); %Open circuit voltages across degradation
PFFs = zeros(num_files,1); %Pseudo fill factors across degradation
PEs = zeros(num_files,1); %Pseudo-Efficiencies across degradation
Ts = zeros(num_files,1); %Temperatures across degradation
Taus = zeros(num_files,1); %Take the maximum of the raw lifetime values for each individual measurement
%Now initialize the raw matrices
seconds_raw = zeros(num_files,125); %For each measurement take raw time measurements
MCDs_raw = zeros(num_files,125); %Take minority carrier density values over time and over degradation
taus_raw = zeros(num_files,125); %Take lifetime values over time and over degradation
%num_files = length(degradation_seconds);
for i=1:num_files
    %Determine target filename
%     filename = strcat(base_sample_name,num2str(degradation_seconds(i)),'s');
%     filename = [base_sample_name num2str(degradation_seconds(i)) 's.xlsm'];
    filename = [base_sample_name num2str(degradation_seconds(i)) '.xlsm'];
%     if i == 1
%         filename = base_sample_name;
%     else
%         filename = strcat(base_sample_name,'_',num2str(degradation_seconds(i)));
%     end
    %Now that we have the target filename and have made the appropriate initializations, use xlsread to take data
    data_summary = xlsread(filename,'Summary','A2:S2');
    data_raw = xlsread(filename,'RawData','A2:J126')'; %Take transpose for syntax simplicity
    %After reading this data from MATLAB, write to vectors for plotting
    PEs(i) = data_summary(2); %efficiency
    VOCs(i) = data_summary(3); %Voc
    PFFs(i) = data_summary(7); %pseudo fill factor 
    Ts(i) = data_summary(17); %temperature
    seconds_raw(i,:) = data_raw(1,:);
    MCDs_raw(i,:) = data_raw(9,:);
    taus_raw(i,:) = data_raw(10,:);
    %We can append to the Taus vector by taking the max of the raw values
    %for this measurement
    Taus(i) = max(data_raw(10,:));
    %Interpolate the MCD tau vectors at the target excess carrier density
    Taus(i) = interp1(MCDs_raw(i,:),taus_raw(i,:),deltan_target); 
end
%After we've taken data, develop plots and files to save
%Plot 1: Injection-Dependent Lifetime
%Now plot excess carrier density 
for i=1:num_files
    loglog(MCDs_raw(i,:),taus_raw(i,:),'LineWidth',1.5);
    hold on
end
title('Injection-Dependent Lifetime \tau','FontSize',15);
set(gca,'FontSize',12);
xlabel('$\displaystyle\frac{C}{m^{3}}$','interpreter','latex','FontSize',15);
ylabel('seconds');
grid on
hold off

%Plot 2: Lifetime over time
figure
plot(degradation_seconds,Taus,'LineWidth',1.5);
set(gca,'FontSize',12);
title('Lifetime \tau vs. t','FontSize',20);
xlabel('seconds','FontSize',15);
ylabel('Lifetime \tau (seconds)', 'FontSize',15);
grid on

%Plot 3: Voc over time
figure
plot(degradation_seconds,VOCs,'LineWidth',1.5);
set(gca,'FontSize',12);
title('V_o_c vs. t','FontSize',20);
xlabel('seconds','FontSize',15);
ylabel('V_o_c (V)', 'FontSize',15);
grid on

%Plot 4: efficiency over time
figure
plot(degradation_seconds,PEs,'LineWidth',1.5);
set(gca,'FontSize',12);
title('Pseudo efficiency \eta vs. t','FontSize',20);
xlabel('seconds','FontSize',15);
ylabel('Pseudo efficiency \eta (%)', 'FontSize',15);
grid on

%Plot 5: temperature over time
figure
plot(degradation_seconds,Ts,'LineWidth',1.5);
set(gca,'FontSize',12);
title('Temperature vs. t','FontSize',20);
xlabel('seconds','FontSize',15);
ylabel('Temperature (C)', 'FontSize',15);
grid on




