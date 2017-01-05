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

%This script analyzes TIDLS measurements taken with WCT-120TS. 
clear all; close all; 
directory = 'C:\Users\Malloryj\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\A41-6\Summary files'; 
before_directory = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160727\for_processing\before';
after_directory = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160727\for_processing\after';
processing_directory = 'C:\Users\Malloryj\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\A41-6\Summary files';
SRV_directory = 'C:\Users\Malloryj\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\16-6-28-P-2\Summary files';
deg_directory = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\PERC LeTID Advanced System\Harmonic difference\degraded';
undeg_directory = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\PERC LeTID Advanced System\Harmonic difference\undegraded';
%type - p or n
type = 'p';
%Fit range for Joe
fit_range = 0.3; 
%% Collect and process raw data - WCT-120TS files
%Find all of the files in the directory
[fileList,fileListShort] = getAllFiles(directory); 
%Get the injection-dependent lifetime data for all files
process_xls_data(directory,[directory '\Raw_data.mat']);
%We also want to store all of the information for each file
%T, thickness, resistivity (entered/measured), type, optical constant, calibration,
%1/64 or 1/1
for file = 1:length(fileList)
    this_file = fileList{file};
    thick{file,1} = xlsread(this_file,'User','B6');
    res{file,1} = xlsread(this_file,'User','C6');
    oc{file,1} = xlsread(this_file,'User','E6');
    temp{file,1} = xlsread(this_file,'Summary','D2');
    meas_res{file,1} = xlsread(this_file,'Summary','Q2');
    calib{file,1} = xlsread(this_file,'Settings','C5');
    doping{file,1} = xlsread(this_file,'Summary','F2');
end
info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
save([directory '\meas_info.mat'],'info'); 
%% Collect and process raw data - WCT-120 OLD files
%Find all of the files in the directory
[fileList,fileListShort] = getAllFiles(directory); 
%Get the injection-dependent lifetime data for all files
process_xls_data(directory,[directory '\Raw_data.mat']);
%We also want to store all of the information for each file
%T, thickness, resistivity (entered/measured), type, optical constant, calibration,
%1/64 or 1/1
for file = 1:length(fileList)
    this_file = fileList{file};
    thick{file,1} = xlsread(this_file,'User','B6');
    res{file,1} = xlsread(this_file,'User','C6');
    oc{file,1} = xlsread(this_file,'User','E6');
    temp{file,1} = 25;
    meas_res{file,1} = xlsread(this_file,'Summary','Q2');
    calib{file,1} = xlsread(this_file,'Summary','T2');
    doping{file,1} = xlsread(this_file,'Summary','E2');
end
info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
save([directory '\meas_info.mat'],'info'); 

%% Collect and process raw data - WCT-120TS temperature sweep
%Note - this script just calculates the data which is stored in the "plots"
%spreadsheet. The actual temperatures are different from what is recorded
%here, and there are more lifetime curves available in the "multi"
%spreadsheet. However, this raw data need to be processed. 
%Find all of the files in the directory
[fileList,fileListShort] = getAllFiles(directory); 
%All of the data is stored in the "Plots" spreadsheet. Let's assume that we
%have one filename. 
all_data = xlsread(fileList{1},'For analysis'); 
%The second row is a header
% all_data(2,:) = []; 
[points,columns] = size(all_data); 
temps = [2:2:(columns)];
for i = 1:length(temps)
    %We want this to be the same structure as our other script
    deltan = all_data(2:end,temps(i)-1); 
    lifetime = all_data(2:end,temps(i)); 
    figure;
    loglog(deltan,lifetime,'.');
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (seconds)','FontSize',20);
    title(num2str(all_data(1,temps(i))),'FontSize',20);
    dataSave{i} = [deltan,lifetime];
    %Save the temperature because that is not common to all samples
    temp{i,1} = all_data(1,temps(i)); 
end
save([directory '\Raw_data.mat'],'fileListShort','dataSave');
%We also want to store all of the information for each file
%T, thickness, resistivity (entered/measured), type, optical constant, calibration,
%1/64 or 1/1
for file = 1:length(fileList)
    this_file = fileList{file};
    thick{file,1} = xlsread(this_file,'User','B6');
    res{file,1} = xlsread(this_file,'User','C6');
    oc{file,1} = xlsread(this_file,'User','E6');
    meas_res{file,1} = xlsread(this_file,'Summary','Q2');
    calib{file,1} = xlsread(this_file,'Settings','C5');
    doping{file,1} = xlsread(this_file,'Summary','F2');
end
info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
save([directory '\meas_info.mat'],'info'); 

%% Now that the raw data has been processed, plot the curves together
load([directory '\Raw_data.mat']);
load([directory '\meas_info.mat']); 
h=figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(dataSave)
    dataNow = dataSave{i}; 
    curves(i)=loglog(dataNow(:,1),dataNow(:,2),'.'); 
    hold all; 
    label{i} = [info(i).filename ', T=' num2str(info(i).temperature)];
end
%% On the same plot as the previous section, plot other room temperature measurements at MIT before T ramp
data= load([before_directory '\Raw_data.mat']); 
dataSave_before = data.dataSave; 
fileListShort_before = data.fileListShort;
for i = 1:length(dataSave_before)
    dataNow = dataSave_before{i}; 
    curves(length(dataSave)+i)=loglog(dataNow(:,1),dataNow(:,2),'k-'); 
    hold all; 
    label{length(dataSave)+i} = [fileListShort_before{i} ', T=' num2str(25)];
end
%% On the same plot as the previous section, plot other room temperature measurements at MIT after T ramp
data= load([after_directory '\Raw_data.mat']); 
dataSave_after = data.dataSave; 
fileListShort_after = data.fileListShort;
for i = 1:length(dataSave_after)
    dataNow = dataSave_after{i}; 
    curves(length(dataSave)+length(dataSave_before)+i)=loglog(dataNow(:,1),dataNow(:,2),'k--'); 
    hold all; 
    label{length(dataSave)+length(dataSave_before)+i} = [fileListShort_after{i} ', T=' num2str(25)];
end
%% Label and save the current summary plot
xlabel('Excess carrier density [cm^-^3]','FontSize',20); 
ylabel('Lifetime [s]','FontSize',20);
legend(curves',label');
hgsave(h,[directory '\Lifetime w T']);
print(h,'-dpng','-r0',[directory '\Lifetime w T.png']); 
%% Load the data that will be processed and calculate the Joe and lifetime contributions
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
temps_Joe = figure;
norm_temps_Joe = figure; 
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %Now that we have the revised data, fit the Joe
    [Joe,Joe_select{i,1}] = fit_Joe_curve(deltan_rev,tau_rev,info(i).doping,info(i).temperature+273.15,type,fit_range,info(i).thickness); 
    %Now that we have the Joe, we can calculate the surface lifetime
    [tau_surf] = calculate_tau_surf_Joe(Joe_select{i,1},type,info(i).doping,deltan_rev,info(i).thickness,info(i).temperature+273.15);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    for j = 1:length(deltan_rev)
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
    end
    %Finally, we can calculate the SRH lifetime!
    tau_SRH = ((1./tau_rev)-(1./tau_intr)-(1./tau_surf)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6);
    hold all; 
    loglog(deltan_rev,tau_surf.*1e6); 
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6);
    hold all;
    loglog(deltan_rev,tau_SRH.*1e6);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Surface','Intrinsic','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_surf_store{i,1} = tau_surf;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
    figure(temps_Joe); 
    plot(info(i).temperature,Joe_select{i,1},'bo');
    hold on; 
    %I also want to plot the normalized quantity Joe/q*ni2 because this
    %should be constant with T
    if info(i).temperature+273.15 <= 303 & info(i).temperature+273.15 >= 297
        Eg = 1.1242; %eV
        [NC,NV] = DOS_std(info(i).temperature+273.15); %cm^-3  
    else
        [Eg] = Sze(info(i).temperature+273.15); %eV
        %Density of states effective mass model for consistency
        [NC,NV] = DOS_em(info(i).temperature+273.15); %cm^-3  
    end
    %Boltzmann constant
    k_B = 8.61733238e-5; %eV/K  
    ni2 = NC*NV*exp(-Eg/(k_B*(info(i).temperature+273.15))); %cm^-6
    q = 1.602e-19; %C
    norm_quant = Joe_select{i,1}/(q*ni2); 
    figure(norm_temps_Joe)
    plot(info(i).temperature,norm_quant,'ro');
    hold on; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_surf',tau_surf_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'Joe',Joe_select);
figure(temps_Joe);
xlabel('Temperature [C]','FontSize',20);
ylabel('J_o_e [A/cm^2]','FontSize',20);
figure(norm_temps_Joe);
xlabel('Temperature [C]','FontSize',20);
ylabel('J_o_e/qn_i^2','FontSize',20);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
%% Load the data that will be processed and calculate lifetime contributions - Joe set by an average
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
cutoff(1,1) = 1e14;
cutoff(1,2) = 1e16;
%fit parameters for Joe as function of T (ax^2+bx+c). T is fit in units of
%C, Joe in units A/cm2. 
a = -2.8082e-4;
b = 1.6869e-1;
c = -3.4422e1;
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
%     disp('Select the region for cutting off the HIGH injection data');
%     [cutoff,nothing]=ginput(1);
%     [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff(1,2));
    %We might always want to remove some low injection data
%     disp('Select the region for cutting off the LOW injection data');
%     [cutoff,nothing]=ginput(1);
%     [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff(1,1));
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %Now that we have the revised data, acquire the Joe from our fitted
    %curve
    ln_Joe = (a)*((info(i).temperature)^2)+((b)*info(i).temperature)+c;  
    Joe{i,1} = exp(ln_Joe); 
    %Now that we have the Joe, we can calculate the surface lifetime
    [tau_surf] = calculate_tau_surf_Joe(Joe{i,1},type,info(i).doping,deltan_rev,info(i).thickness,info(i).temperature+273.15);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    for j = 1:length(deltan_rev)
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
    end
    %Finally, we can calculate the SRH lifetime!
    tau_SRH = ((1./tau_rev)-(1./tau_intr)-(1./tau_surf)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6);
    hold all; 
    loglog(deltan_rev,tau_surf.*1e6); 
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6);
    hold all;
    loglog(deltan_rev,tau_SRH.*1e6);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Surface','Intrinsic','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_surf_store{i,1} = tau_surf;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_surf',tau_surf_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'Joe',Joe);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
%% Load the data that will be processed and calculate lifetime contributions - Joe and cutoffs pre-defined
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
% Joe{1,1} = 5.75e-14.*1;
% % Joe{2,1} = 4.76e-14.*1; 
% cutoff(1,1) = 1e12;
% cutoff(1,2) = 1e18;
% cutoff(2,1) = 1e12;
% cutoff(2,2) = 1e18;
Joe_in = 4.76e-14*0.7; 
cutoff = ones(length(dataSave),2); 
cutoff(:,1) = cutoff(:,1).*2e14; 
cutoff(:,2) = cutoff(:,2).*1e16; 
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    Joe{i,1} = Joe_in; 
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff(i,2));
    %We might always want to remove some low injection data
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff(i,1));
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %Now that we have the Joe, we can calculate the surface lifetime
    [tau_surf] = calculate_tau_surf_Joe(Joe_in,type,info(i).doping,deltan_rev,info(i).thickness,info(i).temperature+273.15);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    for j = 1:length(deltan_rev)
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
    end
    %Finally, we can calculate the SRH lifetime!
    tau_SRH = ((1./tau_rev)-(1./tau_intr)-(1./tau_surf)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6);
    hold all; 
    loglog(deltan_rev,tau_surf.*1e6); 
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6);
    hold all;
    loglog(deltan_rev,tau_SRH.*1e6);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Surface','Intrinsic','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_surf_store{i,1} = tau_surf;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_surf',tau_surf_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'Joe',Joe);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');

%% Load the data that will be processed and lifetime contributions - no surface contribution
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    for j = 1:length(deltan_rev)
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
    end
    %Finally, we can calculate the SRH lifetime, assuming no surface
    %contribution
    tau_SRH = ((1./tau_rev)-(1./tau_intr)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6);
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6);
    hold all;
    loglog(deltan_rev,tau_SRH.*1e6);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Intrinsic','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');

%% Load the data and calculate the SRV as a function of injection for each curve
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']);
SRVtoSave = cell(length(dataSave),1);
h=figure; 
label = zeros(length(dataSave),1)
%Loop through each temperature calculate the SRV
for i = 1:length(dataSave)
    dataNow = dataSave{i}; 
    tau = dataNow(:,2); 
    deltan = dataNow(:,1); 
    %We might want to trim the lifetime
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %We need the intrinsic lifetime
    deltan = deltan_rev;
    tau = tau_rev; 
    tau_intr = zeros(size(deltan));
    diffusivity_save = zeros(size(deltan)); 
    for j = 1:length(deltan)
        %Get the intrinsic lifetime
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan(j),info(i).doping,type);
        %unfortunately the diffusivity is also temperature and injection
        %dependent
        [De,Dh] = diffusivity(info(i).temperature+273.15,type,info(i).doping,deltan(j));
        if type == 'n'
            diffusivity_save(j,1) = Dh; %hole is minority carrier
        elseif type == 'p'
            diffusivity_save(j,1) = De; %electron is minority carrier
        end
    end
    %Calculate the surface-related lifetime assuming zero SRH contribution
    tau_surf = ((1./tau)-(1./tau_intr)).^(-1);
    %Calculate the SRV, including the injection-dependent diffusivity
    SRV = info(i).thickness./((tau_surf-((1./diffusivity_save).*((info(i).thickness/pi)^2))).*2);
    SRVtoSave{i} = [deltan SRV]; 
    %Plot the result so that we can compare
    figure(h);
    loglog(deltan,SRV,'LineWidth',2);
    label(i,1)=info(i).temperature; 
    hold all; 
end
xlabel('Excess carrier density [cm^-^3]','FontSize',20); 
ylabel('Surface recombination velocity [cm/s]','FontSize',20);
legend(num2str(label));
hgsave(h,[directory '\SRV w T']);
print(h,'-dpng','-r0',[directory '\SRV w T.png']); 
save([directory '\SRV_data.mat'],'SRVtoSave'); 

%% Load the data and process the lifetime contributions with SRV given by a .mat file
%Load the data for processing
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
% %Read the temperature dependent diffusivity 
% data = xlsread([directory '\diffusivity.xlsx']); 
%Load the SRV data
load([SRV_directory '\SRV_data.mat']); 
% if type == 'n'
%     %The minority carrier is the hole
%     diffusivity = data(:,5); 
% elseif type == 'p'
%     %the minority carrier is the electron
%     diffusivity = data(:,6); 
% end
% SRVindex = data(:,7); 
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    diffusivity_save = zeros(size(deltan)); 
    for j = 1:length(deltan_rev)
        %Get the intrinsic lifetime
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
        %unfortunately the diffusivity is also temperature and injection
        %dependent
        [De,Dh] = diffusivity(info(i).temperature+273.15,type,info(i).doping,deltan_rev(j));
        if type == 'n'
            diffusivity_save(j,1) = Dh; %hole is minority carrier
        elseif type == 'p'
            diffusivity_save(j,1) = De; %electron is minority carrier
        end
    end
    %We need the surface contribution which comes from a reference wafer
%     SRV_now = SRVtoSave{SRVindex(i)}; 
    SRV_now = SRVtoSave{i};
    %Interpolate the SRV so that it matches the measured lifetime
    SRVq = interp1(SRV_now(:,1),SRV_now(:,2),deltan_rev); 
    tau_surf =(info(i).thickness./(2.*SRVq))+((1./diffusivity_save(i)).*((info(i).thickness/pi)^2)); %cm/s
    %Finally, we can calculate the SRH lifetime
    tau_SRH = ((1./tau_rev)-(1./tau_intr)-(1./tau_surf)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6,'LineWidth',2);
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6,'LineWidth',2);
    hold all;
    loglog(deltan_rev,tau_surf.*1e6,'LineWidth',2); 
    hold all; 
    loglog(deltan_rev,tau_SRH.*1e6,'LineWidth',2);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Intrinsic','Surface','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
    tau_surf_store{i,1} = tau_surf; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'tau_surf',tau_surf_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');

%% Load the data and ask for where to crop the SRH data based on the contributions
%Load the data for processing
load([processing_directory '\lifetime_breakdown.mat']);
for i = 1:length(dataSave);
    deltan_rev = lifetime_breakdown(i).deltan; 
    tau_rev = lifetime_breakdown(i).tau; 
    tau_SRH = lifetime_breakdown(i).tau_SRH; 
    tau_intr = lifetime_breakdown(i).tau_intr;
    tau_surf = lifetime_breakdown(i).tau_surf;
    figure;
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6,'LineWidth',2);
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6,'LineWidth',2);
    hold all;
    loglog(deltan_rev,tau_surf.*1e6,'LineWidth',2); 
    hold all; 
    loglog(deltan_rev,tau_SRH.*1e6,'LineWidth',2);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Intrinsic','Surface','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_highinj(deltan_rev,tau_SRH,cutoff);
    [deltan_rev,tau_rev_rev] = remove_highinj(deltan_rev,tau_rev,cutoff);
    [deltan_rev,tau_surf_rev] = remove_highinj(deltan_rev,tau_surf,cutoff);
    [deltan_rev,tau_intr_rev] = remove_highinj(deltan_rev,tau_intr,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_lowinj(deltan_rev,tau_SRH_rev,cutoff);
    [deltan_rev,tau_rev_rev] = remove_lowinj(deltan_rev,tau_rev_rev,cutoff);
    [deltan_rev,tau_surf_rev] = remove_lowinj(deltan_rev,tau_surf_rev,cutoff);
    [deltan_rev,tau_intr_rev] = remove_lowinj(deltan_rev,tau_intr_rev,cutoff);
     %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev_rev;
    tau_SRH_store{i,1} = tau_SRH_rev;
    tau_intr_store{i,1} = tau_intr_rev; 
    tau_surf_store{i,1} = tau_surf_rev;
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'tau_surf',tau_surf_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
%% Perform the lifetime fitting
%Load the data which has been separated into lifeitme contributions
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']); 
fit_tries = 1e5; 
%This loop may be slow due to fitting
for i = 1:length(lifetime_breakdown)
    %Given the SRH lifetime, we just need to linearize the injection level
    %Plot the Murphy linearization for n-type
    %Get sample parameters at specified temperature
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(info(i).temperature+273.15,info(i).doping,type); 
    %Normalized carrier density
    if type == 'p'
        X = (n0+lifetime_breakdown(i).deltan)./(p0+lifetime_breakdown(i).deltan);
    elseif type == 'n'
        X = (p0+lifetime_breakdown(i).deltan)./(n0+lifetime_breakdown(i).deltan);
    end
    tau_SRH_now = lifetime_breakdown(i).tau_SRH; 
    indices = find(tau_SRH_now<0 | isnan(tau_SRH_now)==1); 
    tau_SRH_now(indices) = [];
    X(indices) = [];
    X_store{i} = X; 
    tau_SRH_store{i} = tau_SRH_now;
    %Write the data to an XLS file for easy copying
    
    xlswrite([processing_directory '\Linearized_data.xlsx'],[X,tau_SRH_now],['Sheet' num2str(i)]); 
    %Now we fit. Fit in microseconds since I'm familiar with good MSE for
    %those units
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_SRH_now.*1e6,info(i).temperature,processing_directory,fit_tries);
    to_write = zeros(6,3); 
    to_write(1:2,1) = one_defect{i,1}';
    twodef = two_defects{i,1}; 
    to_write(1:2,2) = twodef(1,:)';
    to_write(3:4,2) = twodef(2,:)';
    threedef = three_defects{i,1}; 
    to_write(1:2,3) = threedef(1,:)';
    to_write(3:4,3) =threedef(2,:)';
    to_write(5:6,3) = threedef(3,:)';
    xlswrite([processing_directory '\Linearized_data.xlsx'],to_write,['Sheet' num2str(i)],'C1:E6'); 
end
best_fits = struct('one_defect',one_defect,'MSE_one',MSE_one,'two_defects',two_defects,'MSE_two',MSE_two,'three_defects',three_defects,'MSE_three',MSE_three,'all_fits',all_parameters_store,'all_MSE',all_MSE_store);
save([processing_directory '\best_fits.mat'],'best_fits');
save([processing_directory '\original_linearized.mat'],'X_store','tau_SRH_store');
%% Make the E-k curves 
%Given our best fits, now we make some assumptions and create the E-k
%curves at each temperature. 
%We need to load our data first
load([processing_directory '\meas_info.mat']); 
load([processing_directory '\best_fits.mat']);
label = zeros(1,1); 
defect1 = figure;
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
defect2 = figure;
tau_defect1 = figure;
tau_defect2 = figure; 
% defect3 = figure;
%I want these to be in order of increasing temperature
% for i = 1:length(info)
%     T(i) = info(i).temperature;
% end
% %Sort 
% [T,for_iteration] = sort(T); 
for i = 1:length(info)
%     index = for_iteration(i);
    index = i; 
%     best_fit = best_fits(index).one_defect; 
    best_fit = best_fits(index).two_defects;
%     best_fit = best_fits(index).three_defects;
%     Let's sort this to try to match up defects between temperatures
%     [slopes,IX] = sort(best_fit(:,1));
%     best_fit(:,1) = best_fit(IX,1); 
%     best_fit(:,2) = best_fit(IX,2); 
    [m,n] = size(best_fit);
    for j = 1:m
        [Et{index,j},k{index,j},alphanN{index,j}]=generate_Ek(best_fit(j,:),info(index).temperature+273.15,info(index).doping,type);
    end
    figure(defect1); 
    h1(i)=plot(Et{index,1},k{index,1},'-','LineWidth',2,'Color',co{i}); 
    label(i,1) = info(index).temperature; 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{index,1},1./alphanN{index,1},'-','LineWidth',2,'Color',co{i}); 
    hold all;
    figure(defect2);
    h2(i)=plot(Et{index,2},k{index,2},'-','LineWidth',2,'Color',co{i});
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{index,2},1./alphanN{index,2},'-','LineWidth',2,'Color',co{i}); 
    hold all; 
%     figure(defect3);
%     h3(i)=plot(Et{index,3},k{index,3},'-','LineWidth',2,'Color',co{i});
%     hold all;
end
%The Et vector should be the same for each defect
Et_vector = Et{index,j}; 
for i = 1:length(Et_vector)
    [temps,defects] = size(k); 
    for j = 1:defects
        for q = 1:temps
            k_now = k{q,j};
            if isnan(k_now(i))==1
               all_T(j,q) = NaN; 
            else
                all_T(j,q) = k_now(i);
            end
        end
        std_dev(i,j) = std(all_T(j,:)); 
        %We should also track the mean
        average(i,j) = mean(all_T(j,:));
    end
end
figure(defect1); 
%Now we find the minimum of standard deviation curve
% min_std = min(std_dev(:,1)); 
% Et_min = Et_vector(find(std_dev(:,1)==min_std));
% average_k = average(find(std_dev(:,1)==min_std),1); 
% hold all;
% plot([Et_min Et_min], [-100 100],'k--','LineWidth',2);
% hold all;
% plot([0 1.124], [average_k average_k],'k--','LineWidth',2);
% axis([0 1.124 0 100]);
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('Defect 1','FontSize',30); 
figure(defect2); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h2,num2str(label));
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h3,num2str(label));
title('Defect 1','FontSize',30);
figure(tau_defect2); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h4,num2str(label));
title('Defect 2','FontSize',30);
% figure(defect3); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
% ylabel('k [-]','FontSize',20);
% legend(h3,num2str(label));
% title('Defect 3','FontSize',30); 
figure; 
plot(Et_vector,std_dev(:,1),'LineWidth',3);
xlabel('E_t-E_v [eV]','FontSize',20);
ylabel('Standard deviation in k','FontSize',20); 
title('Defect 1','FontSize',30); 
figure; 
plot(Et_vector,std_dev(:,2),'LineWidth',3); 
xlabel('E_t-E_v [eV]','FontSize',20);
ylabel('Standard deviation in k','FontSize',20); 
title('Defect 2','FontSize',30); 
%% Plot all the E-k curves on the same plot because we don't know which defect belongs where
%We need to load our data first
load([processing_directory '\meas_info.mat']); 
load([processing_directory '\best_fits.mat']);
all_defects = figure;
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
tau_defects = figure;
defectcount = 1;
label = zeros(length(best_fits),1);
for i = 1:length(best_fits)
    index = i;
%     best_fit = best_fits(index).one_defect;
    best_fit = best_fits(index).two_defects;
%     best_fit = best_fits(index).three_defects; 
    [m,n] = size(best_fit);
    for j = 1:m
        [Et{index,j},k{index,j},alphanN{index,j}]=generate_Ek(best_fit(j,:),info(index).temperature+273.15,info(index).doping,type);
    end
    figure(all_defects); 
    h1(defectcount)=plot(Et{index,1},k{index,1},'-','LineWidth',2,'Color',co{i}); 
    label(defectcount,1) = info(index).temperature; 
    hold all; 
    figure(tau_defects); 
    h2(defectcount)=plot(Et{index,1},1./alphanN{index,1},'-','LineWidth',2,'Color',co{i}); 
    hold all;
    defectcount = defectcount+1;
    figure(all_defects);
    plot(Et{index,2},k{index,2},'--','LineWidth',2,'Color',co{i});
    hold all;
    figure(tau_defects); 
    plot(Et{index,2},1./alphanN{index,2},'--','LineWidth',2,'Color',co{i}); 
    hold all; 
%     defectcount = defectcount+1;
%     figure(all_defects);
%     plot(Et{index,3},k{index,3},'-.','LineWidth',2,'Color',co{i});
%     hold all;
%     figure(tau_defects); 
%     plot(Et{index,3},1./alphanN{index,3},'-.','LineWidth',2,'Color',co{i}); 
%     defectcount = defectcount+1;
end
figure(all_defects); 
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('All defects','FontSize',30); 
figure(tau_defects); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h2,num2str(label));
title('All defects','FontSize',30);

%% Make E-k curves for room temperature with simple values
%We need to load our data first
load([processing_directory '\meas_info.mat']); 
% load([processing_directory '\best_fits.mat']);
best_fits(2).two_defects = [9.70374e-6, 0.000122852;0.00018, 5.23269e-6];
best_fits(1).two_defects = [0.000150925 8.93748e-5; 0.005911841, 1e-12];
for i = 1:length(info)
    best_fit_now = best_fits(i).two_defects;
    [m,n] = size(best_fit_now);
    for j = 1:m
        best_fit = best_fit_now(j,:); 
        %Boltzmann constant
        k_B = 8.61733238e-5; %eV/K  
        [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(info(i).temperature+273.15,info(i).doping,type); 
        NC = 3e19;
        NV = 1e19;
        Eg = 1.1242;
        vth_e = 2.05e7; 
        vth_h = 1.69e7; 
        T = info(i).temperature+273.15;
        %Define the energy levels for evaluation
        Et = linspace(0,Eg,250); %eV
        Q = zeros(size(Et)); 
        alphanN = zeros(size(Et));
        k = zeros(size(Et)); 
        A = best_fit(1)+best_fit(2); %X -> 1
        B = best_fit(2); %X -> 0
        C = best_fit(1)/A; %slope/X -> 1
        for l = 1:length(Et)
            %Calculate n1
            n1 = NC*exp(-(Eg-Et(l))/(k_B*T)); 
            %Calculate p1
            p1 = NV*exp(-Et(l)/(k_B*T)); 
            %Calculate the Q values for these defects
            Q(l) = (C+(p1/p0))/(1-(n1/p0)-C);
            %Calculate the quantity alphan*Nt for these defects
            alphanN(l) = (1/B)*(1+((1/p0)*((Q(l)*n1)+p1))); 
            %Calculate the k values for these defects
            k(l) = Q(l)*vth_h/vth_e; 
        end
        figure; 
        plot(Et,k); 
        xlabel('E_t-E_v [eV]','FontSize',20); 
        ylabel('k [-]','FontSize',20);
        title([info(i).filename ', Defect ' num2str(j)])
    end
end
%% Perform TDLS analysis based on low injection information
%Load the data which has been separated into lifetime contributions
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']); 
%This loop may be slow due to fitting
for i = 1:length(lifetime_breakdown)
    deltan = lifetime_breakdown(i).deltan;
    tau_SRH = lifetime_breakdown(i).tau_SRH; 
    %Low level injection is at the end of the curve. 
    deltan_LLI = deltan(length(deltan)-5:end); 
    tau_SRH_LLI = tau_SRH(length(deltan)-5:end);
    %Let's plot what we're choosing
    figure;
    loglog(deltan,tau_SRH.*1e6,'.'); 
    hold all;
    loglog(deltan_LLI,tau_SRH_LLI.*1e6,'x'); 
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('SRH','LLI points');
    title(['Temperature = ' num2str(info(i).temperature)]);
    tau_LLI(i) = mean(tau_SRH_LLI); 
    T(i) = info(i).temperature+273.15;
end
%Plot the results
figure;
semilogy(1000./T,tau_LLI.*1e6,'o'); 
xlabel('1000/T [1/K]','FontSize',20);
ylabel('Lifetime in low-level injection (\mus)','FontSize',20);   
%% Make best fits from Excel 
filename = [processing_directory '\IDLS two and three curve fitting.xlsm'];
load([processing_directory '\meas_info.mat']); 
fits = xlsread(filename,'Summary','B3:E11'); %Note - they should be sorted in increasing temperature
[m,n] = size(fits); 
% load([processing_directory '\best_fits.mat']); %We will load this and then replace two defects fits
% for i = 1:length(info)
%     T(i) = info(i).temperature;
% end
% %Sort 
% [T,for_iteration] = sort(T); 
for i = 1:m
%     index = for_iteration(i);
    index = i;
    defect_1 = fits(i,1:2); 
    defect_2 = fits(i,3:4); 
%     best_fits(index).two_defects = [defect_1;defect_2];
    two_defects{i,1} = [defect_1;defect_2];
end
best_fits = struct('two_defects',two_defects);
save([processing_directory '\best_fits.mat'],'best_fits');
%% Plot the lifetime as a function of T which is analyzed
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']); 
h=figure('units','normalized','outerposition',[0 0 1 1]);
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
for i = 1:length(info)
    T(i) = info(i).temperature;
end
%Sort 
[T,for_iteration] = sort(T); 
for i = 1:length(lifetime_breakdown)
    index = for_iteration(i);
    deltan = lifetime_breakdown(index).deltan;
    tau = lifetime_breakdown(index).tau; 
    curves(i)=loglog(deltan,tau.*1e6,'.','MarkerSize',15,'Color',co{i}); 
    hold all; 
    label(i) = round(info(index).temperature);
end
xlabel('Excess carrier density [cm^-^3]','FontSize',30); 
ylabel('Lifetime [\mus]','FontSize',30);
legend(curves,num2str(label));
%% Make the E-k curves for roomT data, varying surface contribution
%Given our best fits, now we make some assumptions and create the E-k
%curves at each temperature. 
%We need to load our data first
load([processing_directory '\meas_info.mat']); 
load([processing_directory '\best_fits.mat']);
defect1 = figure;
defect2 = figure;
tau_defect1 = figure;
tau_defect2 = figure; 
for i = 1:length(best_fits)
    index = i;
    best_fit = best_fits(index).two_defects;
    for j = 1:length(best_fit)
        [Et{index,j},k{index,j},alphanN{index,j}]=generate_Ek(best_fit(j,:),info(1).temperature+273.15,info(1).doping,type);
    end
    figure(defect1);
    if mod(index,2)==0
        %degraded
        co = [0 0 0];
    else
        co = [0 0 1];
    end
    if index >= 3
        line = '--'; 
    else
        line = '-';
    end
    h1(i)=plot(Et{index,1},k{index,1},line,'LineWidth',2,'Color',co); 
    label(i,1) = info(1).temperature; 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{index,1},1./alphanN{index,1},line,'LineWidth',2,'Color',co); 
    hold all;
    figure(defect2);
    h2(i)=plot(Et{index,2},k{index,2},line,'LineWidth',2,'Color',co);
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{index,2},1./alphanN{index,2},line,'LineWidth',2,'Color',co); 
    hold all; 
end
figure(defect1); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('Defect 1','FontSize',30); 
figure(defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h2,num2str(label));
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h3,num2str(label));
title('Defect 1','FontSize',30);
figure(tau_defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h4,num2str(label));
title('Defect 2','FontSize',30);
%% Plot the k-values for Ti-dd for each temperature on top of current plot
%Assume a constant defect energy level
Et_Tidd = 0.28; %eV
for i = 1:length(info)
%     index = for_iteration(i);
    index = i; 
    temp = info(index).temperature;%celcius
    k_Tidd(i) = calculate_Tidd(temp);
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(temp+273.15,info(index).doping,type); 
    figure(defect1);
%     figure(all_defects)
    hold all;
    plot((Et_Tidd-Eiv),k_Tidd(i),'o','MarkerSize',10,'Color',co{i});
%     figure(defect2);
%     hold all;
%     plot(Et_Tidd,k_Tidd(i),'o','MarkerSize',10,'Color',co{i});
end
%% Plot the k-values for Mo-d for each temperature on top of current plot
%Assume a constant defect energy level
Et_Mod_min = 0.28; %eV
Et_Mod_max = 0.375; %eV
for i = 1:length(info)
%     index = for_iteration(i);
    index = i; 
    temp = info(index).temperature;%celcius
    k_Mod(i) = calculate_Mod(temp+273.15);
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(temp+273.15,info(index).doping,type); 
%     figure(all_defects)
    figure(defect1);
    hold all;
    plot([(Et_Mod_min-Eiv) (Et_Mod_max-Eiv)],[k_Mod(i) k_Mod(i)],'x-','MarkerSize',10,'Color',co{i});
%     figure(defect2);
%     hold all;
%     plot([Et_Mod_min Et_Mod_max],[k_Mod(i) k_Mod(i)],'x-','MarkerSize',10,'Color',co{i});
end

%% Try plotting graphs to deduce upper/lower bandgap half per Paudyal
%Load the data which has been separated into lifeitme contributions
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']); 
upper_bgh = figure;
eval_deltan=1e15;
lower_bgh = figure; 
for i = 1:length(lifetime_breakdown)
    deltan = lifetime_breakdown(i).deltan;
    tau_SRH = lifetime_breakdown(i).tau_SRH; 
    eval_tau_SRH = interp1(deltan,tau_SRH,eval_deltan); 
    k_B = 8.61733238e-5; %eV/K  
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(info(i).temperature+273.15,info(i).doping,type); 
    [NC,NV] = DOS_em(info(i).temperature+273.15); %cm^-3
    [Eg] = Sze(info(i).temperature+273.15); %eV
    [vth_e,vth_h] = vth_em(info(i).temperature+273.15); %cm/s
    %Titanium hypothesis
    n1 = NC*exp(-0.27/(k_B*(info(i).temperature+273.15)));
    x_value = (n1+eval_deltan)./(p0+eval_deltan);
    figure(upper_bgh);
    plot(x_value,eval_tau_SRH,'b.','MarkerSize',10); 
    hold on; 
    deltan_lowinj = find(deltan<=3e14); 
    tau_SRH_lowinj = tau_SRH(deltan_lowinj); 
    deltan_lowinj = deltan(deltan_lowinj); 
    figure(lower_bgh); 
    hold all;
    loglog(deltan_lowinj,tau_SRH_lowinj,'-o','LineWidth',3'); 
    labels(i) = info(i).temperature;
end
figure(upper_bgh); 
xlabel('(n_1(T)+\Deltan)/(p_0(t)+\Deltan)','FontSize',20);
ylabel('\tau_{SRH}','FontSize',20); 
title(['\Deltan = ' num2str(eval_deltan)],'FontSize',20); 
figure(lower_bgh); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('\tau_{SRH}','FontSize',20); 
legend(num2str(labels'));

%% Make SRH lifetime from harmonic sum
%Load the data
load([undeg_directory '\Raw_data.mat']);
data_undeg = dataSave;
load([undeg_directory '\meas_info.mat']); 
info_undeg = info; 
load([deg_directory '\Raw_data.mat']);
data_deg = dataSave;
load([deg_directory '\meas_info.mat']); 
info_deg = info; 
SRH_fig = figure; 
labels = cell(length(data_undeg),1); 
%We assume that they are in the same order in the .mat files
for i = 1:length(data_undeg)
    %Plot the two lifetimes together. 
    h=figure; 
    dataNow_undeg = data_undeg{i};
    dataNow_deg = data_deg{i}; 
    loglog(dataNow_undeg(:,1),dataNow_undeg(:,2),'.'); 
    hold all; 
    loglog(dataNow_deg(:,1),dataNow_deg(:,2),'.'); 
    xlabel('Excess carrier density [cm^-^3]','FontSize',30); 
    ylabel('Lifetime [s]','FontSize',30);
%     legend('Undegraded','Degraded');
    legend('PDG','As-grown');
    title(['Temperature = ' num2str(info_undeg(i).temperature) 'C, ' num2str(info_deg(i).temperature) 'C'],'FontSize',30);
    hgsave(h,[directory '\Lifetime' num2str(round(info_undeg(i).temperature))]);
    print(h,'-dpng','-r0',[directory '\Lifetime ' num2str(round(info_undeg(i).temperature)) '.png']); 
    %Now calculate the SRH lifetime from the two contributions
    %First we need to interpolate. 
    %But before we do that we need to make sure there aren't any repeated
    %injection levels
    [deltan_undeg,tau_undeg] = remove_duplicates(dataNow_undeg(:,1),dataNow_undeg(:,2));
    [deltan_deg,tau_deg] = remove_duplicates(dataNow_deg(:,1),dataNow_deg(:,2));
    tau_undeg_interp = interp1(deltan_undeg,tau_undeg,deltan_deg);
    %Quick correction in case the "degraded" is lower lifetime than the
    %undegraded
    if interp1(deltan_undeg,tau_undeg,1e15)>interp1(deltan_deg,tau_deg,1e15)
        %This is the situation we expect
        %Now everything should be at the same injection levels
        tau_SRH_now = ((1./tau_deg)-(1./tau_undeg_interp)).^(-1);
    elseif interp1(deltan_undeg,tau_undeg,1e15)<interp1(deltan_deg,tau_deg,1e15)
        disp(['Degraded lifetime is higher than undegraded at T = '...
            num2str(info_undeg(i).temperature) ', ' num2str(info_deg(i).temperature)]);
        tau_SRH_now = ((1./tau_undeg_interp)-(1./tau_deg)).^(-1);
    end
    %Plot the result along with the other temperatures
    figure(SRH_fig); 
    curves(i) = loglog(deltan_deg,tau_SRH_now,'LineWidth',2);
    hold all; 
    labels{i} = [num2str(round(info_undeg(i).temperature))];
    %we want to put this in the format which is expected
    tau_SRH_store{i,1} = tau_SRH_now; 
    deltan_store{i,1} = deltan_deg; 
end
%Label the figure
xlabel('Excess carrier density [cm^-^3]','FontSize',30); 
ylabel('Lifetime [s]','FontSize',30);
legend(curves',labels);
hgsave(h,[directory '\SRH lifetimes']);
print(h,'-dpng','-r0',[directory '\SRH lifetimes.png']); 
lifetime_breakdown = struct('deltan',deltan_store,'tau_SRH',tau_SRH_store);
% save([deg_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
%% %% Load the data and ask for where to crop the SRH data based on the contributions
%Load the data for processing
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']);
for i = 1:length(lifetime_breakdown);
    deltan_rev = lifetime_breakdown(i).deltan; 
    tau_SRH = lifetime_breakdown(i).tau_SRH; 
    figure;
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_SRH.*1e6,'LineWidth',2);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    title(['Temperature = ' num2str(info(i).temperature)]);
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_highinj(deltan_rev,tau_SRH,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_lowinj(deltan_rev,tau_SRH_rev,cutoff);
     %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_SRH_store{i,1} = tau_SRH_rev;
end
lifetime_breakdown = struct('deltan',deltan_store,'tau_SRH',tau_SRH_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
 
    
