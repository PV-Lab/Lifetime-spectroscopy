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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017'; 
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','C-1','C-2','H-1','H-2','FZ'};
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12','68-2'};
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12','68-2','66-2','68-4'};
% samples = {'66-2','68-4','S1-lg'};
samples = {'FZ','68-2','66-2','68-4'};
% samples = {'56a','60a','61a','H-1','H-2','FZ','FZ-12','68-2'};
% samples = {'FZ-12','68-2'};
% samples = {'H-1','H-2','FZ','FZ-12','68-2'};
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

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017'; 
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','C-1','C-2','H-1','H-2','FZ'};
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12','68-2','66-2','68-4'};
samples = {'FZ','68-2','66-2','68-4'};
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation';
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 9 2017';
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 13 2017'; 
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 17 2017'; 
dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 19 2017';
dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 23 2017';
dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 25 2017';
dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 27 2017';
dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 31 2017';
dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 2 2017';
dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 6 2017';
dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 8 2017';
dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 13 2017';
dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 15 2017';
dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 17 2017';
% dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 1 2017'; 
dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 1 2017\revised calibration';
dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 8 2017';
dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 11 2017';
dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 15 2017';
% dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 17 2017';
dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 21 2017';
dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 23 2017';
dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 30 2017';
dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 3 2017';
dirname23 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 11 2017';
dirname24 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 13 2017';
dirname25 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 20 2017';
dirname26 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 21 2017';
dirname27 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 27 2017';
dirname28 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 3 2017';
dirname29 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 5 2017';
dirname30 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 9 2017';
dirname31 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 12 2017';
dirname32 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 19 2017';
dirname33 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 25 2017';
dirname34 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 1 2017';
dirname35 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 7 2017';
dirname36 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 13 2017';
dirname37 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\July 12 2017';
dirname38 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017';
% dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 dirname8 dirname9 dirname10 dirname11 dirname12}; 
dirnames = {dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 dirname8 ...
    dirname10 dirname11 dirname12 dirname13 dirname14 dirname15 dirname16 ...
    dirname17 dirname18 dirname19 dirname20 dirname21 dirname22 dirname23 ...
    dirname24 dirname25 dirname26 dirname27 dirname28 dirname29 dirname30 ...
    dirname31 dirname32 dirname33 dirname34,dirname35,dirname36 dirname37 ...
    dirname38}; 
labels = {'initial','1000s','2000s','3000s','4000s','5000s','10000s',...
    '20000s','30000s','40000s','50000s','60000s','70000s','80000s','90000s',...
    '100000s','154495','206005','258325','300025','349990','402730s',...
    '454750s','508360s' '601540s', '702100s','801610s','902260s',...
    '1004410s','1167760s','1408600s','1579630s','1804090s','2024170s',...
    '2509510s','2536210s'};
cm = colormap(hsv(length(dirnames))); 
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','C-1','C-2','H-1','H-2','FZ'};
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12','68-2','66-2','68-4'};
samples = {'FZ','68-2','66-2','68-4'};
savename = '_2536210s_lifetime summary';
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
                t = 2; 
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

%% Make the degradation curves
clear all; close all; clc; 
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation';
savename = '_2536210s_degradation_log';
max_time = 2536210; 
time_shift_E = 801610; %amount of time to shift company E measurements over for comparison after switch
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_removingInitial.xlsx'; 
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

%Now try to get the degradation curve, removing the surface component
SRV_control = 'FZ';
FZ_dop = 5.7e15;
[De,Dh] = diffusivity(300,'p',FZ_dop,deltan_target); 
D_FZ  = De; 
W_FZ = .027; 
index = find(strcmp(SRV_control,samples)==1);
raw_now = lifetime_all{index};
SRV = zeros(length(raw_now),1); 
for i = 1:length(SRV)
    tau_intr = Richter(300,deltan_target,FZ_dop,'p');
    tau_surf = ((1./raw_now(i,2))-(1./tau_intr)).^(-1);
    SRV(i) = W_FZ./((tau_surf-((1/D_FZ)*((W_FZ/pi)^2))).*2);
end
figure;
plot(raw_now(:,1),SRV); 
axis([0 max_time 0 10]);
xlabel('time [s]','FontSize',25); 
ylabel('SRV [cm/s]','FontSize',25); 
SRV_t = raw_now(:,1); 

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
    if raw_now(1,1)==0
        raw_now(1,1) = 1; 
    end
    if norm_now(1,1)==0
        norm_now(1,1) = 1; 
    end
    figure(lifetime_raw); 
    semilogx(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
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
    semilogx(norm_now(:,1),tau_rev./tau_rev(1),'-o','LineWidth',3,'MarkerSize',10);
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
    if raw_now(1,1)==0
        raw_now(1,1) = 1; 
    end
    if norm_now(1,1)==0
        norm_now(1,1) = 1; 
    end
    figure(lifetime_raw); 
    semilogx(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
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
    if raw_now(1,1)==0
        raw_now(1,1) = 1; 
    end
    if norm_now(1,1)==0
        norm_now(1,1) = 1; 
    end
    figure(lifetime_raw); 
    semilogx(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
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
    if raw_now(1,1)==0
        raw_now(1,1) = 1; 
    end
    if norm_now(1,1)==0
        norm_now(1,1) = 1; 
    end
    if strcmp(compE{1,i},'68-2')==1 
        remove = find(raw_now(:,1)<time_shift_E); 
        raw_now(remove,:) = []; 
        remove = find(norm_now(:,1)<time_shift_E); 
        norm_now(remove,:) = []; 
        %Modify the norm calculation for this sample 
        norm_now(:,2) = raw_now(:,2)./raw_now(1,2); 
    end
    raw_now(:,1) = raw_now(:,1)-time_shift_E; 
    norm_now(:,1) = norm_now(:,1)-time_shift_E;
    figure(lifetime_raw); 
    semilogx(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
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

%% A few helpful plots for analysis
% First plot: initial lifetimes versus MIRHP time
% Second plot: Histogram of measured control lifetimes for mc-Si sample,
% with mean and standard deviation
% Third plot: Histogram of measured control lifetimes for FZ passivation 
% sample with mean and standard deviation
% Fourth plot: Histogram of measured control lifetimes for FZ degradation
% sample with mean and standard deviation
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation';
unfired_names = {'61a','54a','50a','45a'};
unfired_H = [0 10 30 120]; 
fired_names = {'49a','53a','56a','52a'};
fired_H = [0 10 30 120]; 
mcSi_controls = {'44a','55a','60a'}; 
mcSi_controls_H = [30 30 0]; 
mcSi_controls_color = {'r','b','k'}; 
Cz_controls = {'H-1','H-2'}; 
Cz_controls_H = [120 120]; 
Cz_controls_color = {'r','b'}; 

start_tau = figure('units','normalized','outerposition',[0 0 1 1]);; 
for i = 1:length(unfired_names)
    index = find(strcmp(unfired_names{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    figure(start_tau); 
    h(1)=plot(unfired_H(i),raw_now(1,2).*1e6,'ro','MarkerSize',10,'LineWidth',2); 
    hold all; 
end
for i = 1:length(fired_names)
    index = find(strcmp(fired_names{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    figure(start_tau); 
    h(2)=plot(fired_H(i),raw_now(1,2).*1e6,'bs','MarkerSize',10,'LineWidth',2); 
    hold all; 
end
for i = 1:length(mcSi_controls)
    index = find(strcmp(mcSi_controls{1,i},samples)==1);
    raw_now = lifetime_all{index}; 
    figure(start_tau); 
    if i == 3
        marker = [mcSi_controls_color{i} 'd'];
        h(4)=plot(mcSi_controls_H(i),raw_now(1,2).*1e6,marker,'MarkerSize',10,'LineWidth',2);
    else
        marker = [mcSi_controls_color{i} '+'];
        h(3)=plot(mcSi_controls_H(i),raw_now(1,2).*1e6,marker,'MarkerSize',10,'LineWidth',2); 
    end
    hold all; 
end
% for i = 1:length(Cz_controls)
%     index = find(strcmp(Cz_controls{1,i},samples)==1);
%     raw_now = lifetime_all{index}; 
%     figure(start_tau); 
%     marker = [Cz_controls_color{i} 'd'];
%     h(5)=plot(Cz_controls_H(i),raw_now(1,2).*1e6,marker,'MarkerSize',10,'LineWidth',2);
% end
figure(start_tau); 
xlabel('MIRHP time [min]','FontSize',25); 
ylabel('initial lifetime [\mus]','FontSize',25); 
% legend(h,'unfired mc-Si','fired mc-Si','mc-Si no H, T only','mc-Si with SiN_x','Cz');
axis([-10 130 80 160]); 
legend(h,'unfired mc-Si','fired mc-Si','mc-Si no H, T only','mc-Si with SiN_x');
title('lifetimes as-received, before degradation','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(start_tau,[savedirname '\initial_lifetimes']);
print(start_tau,'-dpng','-r0',[savedirname '\initial_lifetimes.png']);

%Second plot
mcSi_tau_controls = {'68-2','66-2'}; 
end_times = [801610, max_time]; 
mcSi_controls = figure('units','normalized','outerposition',[0 0 1 1]);; 
for i = 1:length(mcSi_tau_controls)
    index = find(strcmp(mcSi_tau_controls{1,i},samples)==1);
    norm_now = norm_lifetime_all{index};
    indices = find(norm_now(:,1)<= end_times(i)); 
    h(i)=histogram(norm_now(indices,2),10); 
    disp(['Mean for sample ' mcSi_tau_controls{i} ' is ' num2str(mean(norm_now(indices,2)))]); 
    disp(['Standard deviation for sample ' mcSi_tau_controls{i} ' is ' num2str(std(norm_now(indices,2)))]); 
    hold on; 
end
figure(mcSi_controls); 
xlabel('normalized lifetime','FontSize',25); 
ylabel('frequency','FontSize',25); 
title('estimating measurement-to-measurement error','FontSize',25)
legend(h,'68-2','66-2'); 
set(0,'defaultAxesFontSize', 20)
hgsave(mcSi_controls,[savedirname '\mcSi_histogram']);
print(mcSi_controls,'-dpng','-r0',[savedirname '\mcSi_histogram.png']);

%Third plot
FZ_pass_control = {'FZ'}; 
FZ_pass_control_fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(FZ_pass_control)
    index = find(strcmp(FZ_pass_control{1,i},samples)==1);
    norm_now = norm_lifetime_all{index};
    h(i)=histogram(norm_now(:,2),10); 
    disp(['Mean for sample ' FZ_pass_control{i} ' is ' num2str(mean(norm_now(:,2)))]); 
    disp(['Standard deviation for sample ' FZ_pass_control{i} ' is ' num2str(std(norm_now(:,2)))]); 
    hold on; 
end
figure(FZ_pass_control_fig); 
xlabel('normalized lifetime','FontSize',25); 
ylabel('frequency','FontSize',25); 
title('estimating variation in surface passivation quality','FontSize',25)
set(0,'defaultAxesFontSize', 20)
hgsave(FZ_pass_control_fig,[savedirname '\FZ_histogram']);
print(FZ_pass_control_fig,'-dpng','-r0',[savedirname '\FZ_histogram.png']);

%Fourth plot
FZ_deg_control = {'FZ-12'}; 
FZ_deg_control_fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(FZ_deg_control)
    index = find(strcmp(FZ_deg_control{1,i},samples)==1);
    norm_now = norm_lifetime_all{index};
    h(i)=histogram(norm_now(:,2),10); 
    disp(['Mean for sample ' FZ_deg_control{i} ' is ' num2str(mean(norm_now(:,2)))]); 
    disp(['Standard deviation for sample ' FZ_deg_control{i} ' is ' num2str(std(norm_now(:,2)))]); 
    hold on; 
end
figure(FZ_deg_control_fig); 
xlabel('normalized lifetime','FontSize',25); 
ylabel('frequency','FontSize',25); 
title('estimating influence of degradation setup on measurement','FontSize',25)
set(0,'defaultAxesFontSize', 20)
hgsave(FZ_deg_control_fig,[savedirname '\FZ-12_histogram']);
print(FZ_deg_control_fig,'-dpng','-r0',[savedirname '\FZ-12_histogram.png']);