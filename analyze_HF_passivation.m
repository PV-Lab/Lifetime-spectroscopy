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
%% First process the raw data
clear all; close all; clc; 
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\January 03 2018'; 
samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6','P-1'};
for index = 1:length(samples)
    [fileList,fileListShort] = getAllFiles([dirname '\' samples{index}]); 
    savename = [dirname '\' samples{index} '\Raw_data.mat']';
    process_xls_data([dirname '\' samples{index}],savename);
    %We also want to store all of the information for each file
    %T, thickness, resistivity (entered/measured), type, optical constant, calibration,
    %1/64 or 1/1
    thick = cell(size(fileList));
    res = cell(size(fileList));
    oc = cell(size(fileList));
    temp = cell(size(fileList));
    meas_res = cell(size(fileList));
    calib = cell(size(fileList));
    doping = cell(size(fileList));
    for file = 1:length(fileList)
        this_file = fileList{file};
        thick{file,1} = xlsread(this_file,'User','B6');
        res{file,1} = xlsread(this_file,'User','C6');
        oc{file,1} = xlsread(this_file,'User','E6');
        temp{file,1} = 25;
        meas_res{file,1} = xlsread(this_file,'Summary','M2');%'Q2');
        calib{file,1} = xlsread(this_file,'Summary','S2');%'T2');
        doping{file,1} = xlsread(this_file,'Summary','E2');
    end
    info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    save([dirname '\' samples{index} '\meas_info.mat'],'info');
end

%% Now analyze the data
clear all; close all; clc;
%Process data after HF passivation

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\January 03 2018'; 
samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6','P-1'};
lifetime_store = zeros(length(samples),1); 

for i = 1:length(samples)
    load([dirname '\' samples{i} '\Raw_data.mat']); 
    load([dirname '\' samples{i} '\meas_info.mat']); 
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    label = {};
    for j = 1:length(dataSave)
        datanow = dataSave{j}; 
        curves(j)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2); 
        hold all; 
        label{j} = info(j).filename;
        xlim([5e13 1e17])
    end
    xlabel('excess carrier density [cm^-^3]','FontSize',30); 
    ylabel('lifetime [s]','FontSize',30);
    legend(curves',label');
    set(0,'defaultAxesFontSize', 20)
    hgsave(h,[dirname '\' samples{i} '\lifetime summary']);
    print(h,'-dpng','-r0',[dirname '\' samples{i} '\lifetime summary.png']);
    
    %Let's choose measurement 2 and get the lifetime value at 1e15 just for
    %reference
    if length(dataSave)>1
        index = 2; 
    else
        index = 1;
    end
    datanow = dataSave{index}; 
    [deltan,tau] = remove_duplicates(datanow(:,1),datanow(:,2));
    %There's a little bug in this program, for now just do this... 
    try
        lifetime_store(i) = interp1(deltan,tau,1e15); 
    catch
        [deltan,tau] = remove_duplicates(deltan,tau);
        try 
            lifetime_store(i) = interp1(deltan,tau,1e15);
        catch
            [deltan,tau] = remove_duplicates(deltan,tau);
            lifetime_store(i) = interp1(deltan,tau,1e15);
        end
    end
end

%% Analyze different states together
clear all; close all; clc;
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\Summary\50000s';
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 6 2017';
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 8 2017\before';
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 18 2017';
dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 19 2017';
dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\January 03 2018';
dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5}; 
labels = {'initial','1000s','5000s','10000s','50000s'};
cm = colormap(hsv(length(dirnames))); 
samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6','P-1'};
savename = '_50000s_lifetime summary';
for i = 1:length(samples)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    curves = [];
    label = {};
    count = 1;
    for k = 1:length(dirnames)
        try 
            load([dirnames{k} '\' samples{i} '\Raw_data.mat']); 
            load([dirnames{k} '\' samples{i} '\meas_info.mat']); 
            flag = 1; 
        catch
            %If that's not successful, just tell me that the sample doesn't
            %exist for that dirname
            warning(['Error accessing data for directory ' num2str(k) ', sample ' samples{i}]);
            flag = 0; 
        end
        if flag == 1
%             for j = 1:length(dataSave)
            if length(dataSave)>1
                if length(dataSave) == 3
                    t = 2; 
                else
                    t = 1;
                end
            else
                t = 1;
            end
            for j  = t %just the second measurement in each set
                datanow = dataSave{j}; 
                curves(count)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2,'color',cm(k,:)); 
                hold all; 
%                 label{count} = ['Set ' num2str(k) ', #' num2str(j)];
                label{count} = labels{k};
                xlim([5e13 1e17])
                count = count+1; 
            end
        end
    end
    xlabel('excess carrier density [cm^-^3]','FontSize',30); 
    ylabel('lifetime [s]','FontSize',30);
    title(samples{i},'FontSize',30);
    legend(curves',label');
    if strcmp(samples{i},'C-1')==1 || strcmp(samples{i},'C-2')==1 || strcmp(samples{i},'H-1')==1 ||...
        strcmp(samples{i},'H-2')==1 || strcmp(samples{i},'FZ')==1 || strcmp(samples{i},'FZ-12')==1
        ylim([1e-4 1e-2]);
    end
    set(0,'defaultAxesFontSize', 20)
    hgsave(h,[savedirname '\' samples{i} savename]);
    print(h,'-dpng','-r0',[savedirname '\' samples{i} savename '.png']);
end

%% Plot lifetime as function of LTA temperature
clear all; close all; clc; 
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\from SERIS\Lifetime data\summary\initial';
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\from SERIS\Lifetime data\December 5 2017';
samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6'};
LTA = [0 600 550 500 450 700 750 650];
deltan_target = 1e15; 
lifetime_store = zeros(size(LTA)); 

for i = 1:length(samples)
    load([dirname '\' samples{i} '\Raw_data.mat']); 
    load([dirname '\' samples{i} '\meas_info.mat']); 
    if length(dataSave)>1
        if length(dataSave) == 3
            t = 2; 
        else
            t = 1;
        end
    else
        t = 1;
    end
    datanow = dataSave{t}; 
    [deltan,tau] = remove_duplicates(datanow(:,1),datanow(:,2));
    %There's a little bug in this program, for now just do this... 
    try
        lifetime_store(i) = interp1(deltan,tau,deltan_target); 
    catch
        [deltan,tau] = remove_duplicates(deltan,tau);
        try
            lifetime_store(i) = interp1(deltan,tau,deltan_target);
        catch
            [deltan,tau] = remove_duplicates(deltan,tau);
            lifetime_store(i) = interp1(deltan,tau,deltan_target);
        end
    end
    if isnan(lifetime_store(i))==1
        disp(['Lifetime at ' num2str(deltan_target) ' was NaN for sample ' samples{i}]);
    end
end

h=figure;
plot(LTA,lifetime_store.*1e6,'o','MarkerSize',15,'LineWidth',2); 
xlabel('low temperature anneal [C]','FontSize',20); 
ylabel(['lifetime at \Deltan = ' num2str(deltan_target,'%1.0e') ' cm^{-3}'],'FontSize',20);
title('before degradation','FontSize',20); 
set(0,'defaultAxesFontSize', 20)
hgsave(h,[savedirname '\initial lifetime summary with LTA']);
print(h,'-dpng','-r0',[savedirname '\initial lifetime summary with LTA.png']);    

%% Make the degradation curves
clear all; close all; clc; 
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\Summary\50000s';
savename = '_50000s_degradation';
max_time = 50000; 
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\measurement_details.xlsx'; 
deltan_target = 6e14; %target injection level for the measurements, changed to 6e14 on 2/13/17 from 5e14
%Get the measurement details
[meas,samples] = xlsread(meas_details,'measurements');
samples(1,:) = []; 
[times,filenames] = xlsread(meas_details,'filenames'); 
%Make these the same size
filenames = filenames(2:end,2); 

lifetime_all = cell(length(samples),1); 
norm_lifetime_all = cell(length(samples),1); 

for i = 1:length(samples)
    meas_thissample = meas(i,:);
    lifetime_store = [];
    for j = 1:length(meas_thissample)
        if isnan(meas_thissample(j))==0
            %now create the proper filename
            findex = find(meas_thissample(j)==times);  
%             if meas_thissample(j)==0 && length(find(meas_thissample==0))==1
%                 findex=2; %the second round of initial measurements
%             elseif length(find(meas_thissample==0))>1 && meas_thissample(j)==0
%                 findex=findex(j);
%             end
            filename = [filenames{findex} '\' samples{i} '\Raw_data.mat'];
            load(filename);
            if length(dataSave)>1
                t = 2; 
            else
                t = 1;
            end
            datanow = dataSave{t}; 
            [deltan,tau] = remove_duplicates(datanow(:,1),datanow(:,2));
            %There's a little bug in this program, for now just do this... 
            try
                lifetime_store(j,1) = interp1(deltan,tau,deltan_target); 
            catch
                [deltan,tau] = remove_duplicates(deltan,tau);
                try
                    lifetime_store(j,1) = interp1(deltan,tau,deltan_target);
                catch
                    [deltan,tau] = remove_duplicates(deltan,tau);
                    lifetime_store(j,1) = interp1(deltan,tau,deltan_target);
                end
            end
            if isnan(lifetime_store(j,1))==1
                disp(['Lifetime at ' num2str(deltan_target) ' was NaN for sample ' samples{i} ', time ' num2str(meas_thissample(j)) 's']);
            end
        end
    end
    nan_indices = find(isnan(meas_thissample)==1); 
    meas_thissample(nan_indices) = []; 
    lifetime_all{i} = [meas_thissample' lifetime_store]; 
    norm_lifetime_all{i} = [meas_thissample' lifetime_store./lifetime_store(1)];
end

%Now, which samples do we want to plot together?
control = {'H-1','H-2','FZ','FZ-12','68-2','66-2';...
    'Unfired Cz (120 min H)','Fired Cz (120 min H)','FZ passivation','FZ degradation','mc-Si 950C fired undegraded','mc-Si 750C fired undegraded'};
control_param = [5e15 5e15 5e15 5e15 9e15 9e15; .025 .025 .025 .02 .017 0.017]; %first row doping, second row thickness 
fired = {'49a','53a','56a','52a','55a','60a';...
    '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
fired_param = []; %first row doping, second row thickness
unfired = {'61a','54a','50a','45a','44a','60a';...
    '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
unfired_param = []; %first row doping, second row thickness
compE = {'68-2','68-4','66-2';...
    'mc-Si 950C firing no layers','mc-Si 950C firing w/ layers','mc-Si 750C firing undegraded'};
compE_param = []; 

lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm_corr=figure('units','normalized','outerposition',[0 0 1 1]);
[nothing,samp] = size(control); 
labels = {}; 
for i = 1:samp
    index = find(strcmp(control{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_raw); 
    plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
    hold all;  
    %Get the lifetime, accounting for surface
    doping_now = control_param(1,i); 
    W_now = control_param(2,i); 
    [D_now,Dh] = diffusivity(300,'p',doping_now,deltan_target);
    tau_rev = zeros(length(raw_now),1); 
    [num_meas,columns] = size(raw_now); 
    for j = 1:num_meas
        try
            t_index = find(SRV_t==raw_now(j,1)); 
            tau = raw_now(j,2); 
            tau_surf = (W_now./(2.*SRV(t_index)))+((1/D_now).*((W_now/pi)^2)); %cm/s
            tau_rev(j) = ((1./tau)-(1./tau_surf))^(-1); 
        catch
            disp(['There was an error calculating the surface lifetime for ' control{1,i} ', time ' num2str(raw_now(j,1)) 's']);
        end
    end
    figure(lifetime_norm_corr); 
    plot(norm_now(:,1),tau_rev./tau_rev(1),'-o','LineWidth',3,'MarkerSize',10);
    hold all; 
    labels{i,1} = control{2,i};
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
title('control samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\controls' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\controls' savename '.png']);
figure(lifetime_norm); 
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
axis([0 max_time 0 2]);
title('control samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\controls_norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\controls_norm' savename '.png']);
figure(lifetime_norm_corr); 
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime (no surface) [-]','FontSize',25); 
legend(labels); 
axis([0 max_time 0 2]);
title('control samples (surface removed)','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm_corr,[savedirname '\controls_norm_nosurf' savename]);
print(lifetime_norm_corr,'-dpng','-r0',[savedirname '\controls_norm_nosurf' savename '.png']);
    
lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
[nothing,samp] = size(unfired); 
labels = {};
for i = 1:samp
    index = find(strcmp(unfired{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_raw); 
    plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
    hold all; 
    labels{i,1} = unfired{2,i}; 
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
title('unfired samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\unfired' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\unfired' savename '.png']);
figure(lifetime_norm); 
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
axis([0 max_time 0 2]);
title('unfired samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\unfired_norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\unfired_norm' savename '.png']);

lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
[nothing,samp] = size(fired); 
labels = {};
for i = 1:samp
    index = find(strcmp(fired{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_raw); 
    plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
    hold all; 
    labels{i,1} = fired{2,i}; 
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
title('fired samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\fired' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\fired' savename '.png']);
figure(lifetime_norm); 
axis([0 max_time 0 2]);
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
title('fired samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\fired_norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\fired_norm' savename '.png']);

lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
[nothing,samp] = size(compE); 
labels = {};
for i = 1:samp
    index = find(strcmp(compE{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    norm_now = norm_lifetime_all{index};
    if strcmp(compE{1,i},'68-2')==1 
        remove = find(raw_now(:,1)<time_shift_E); 
        raw_now(remove,:) = []; 
        remove = find(norm_now(:,1)<time_shift_E); 
        norm_now(remove,:) = []; 
    end
    raw_now(:,1) = raw_now(:,1)-time_shift_E; 
    norm_now(:,1) = norm_now(:,1)-time_shift_E; 
    figure(lifetime_raw); 
    plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
    hold all; 
    labels{i,1} = compE{2,i}; 
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
title('Company E samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\compE' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\compE' savename '.png']);
figure(lifetime_norm); 
axis([0 max_time 0 2]);
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
title('Company E samples','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\compE_norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\compE_norm' savename '.png']);