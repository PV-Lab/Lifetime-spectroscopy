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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 28 2017'; 
samples = {'44b','45b','49b','50b','52b','53b','54b','55b','56b','60b','61b','C-1','C-2','66-2','68-2','68-4','FZ-12','FZ-new'};
% samples = {'FZ-12'};
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
        meas_res{file,1} = xlsread(this_file,'Summary','N2');%'Q2');
        calib{file,1} = xlsread(this_file,'Summary','T2');%'T2');
        doping{file,1} = xlsread(this_file,'Summary','E2');
    end
    info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    save([dirname '\' samples{index} '\meas_info.mat'],'info');
end

%% Now analyze the data
clear all; close all; clc;
%Process data after HF passivation

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 28 2017'; 
samples = {'44b','45b','49b','50b','52b','53b','54b','55b','56b','60b','61b','C-1','C-2','66-2','68-2','68-4','FZ-12','FZ-new'};
% samples = {'FZ-12'};
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
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017'; 
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 11 2017';
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 16 2017';
dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 29 2017';
dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 30 2017';
dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 8 2017';
dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 14 2017';
dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 20 2017';
dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 22 2017';
dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 28 2017';
samples = {'44b','45b','49b','50b','52b','53b','54b','55b','56b','60b','61b','C-1','C-2','FZ','66-2','68-2','68-4','FZ-12','FZ-new'};
dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
    dirname8 dirname9 dirname10}; 
labels = {'initial' '50s' '100s' '250s' '500s' '501s' '750s' '1000s' '1500s' ...
    '2000s'};
cm = colormap(hsv(length(dirnames)));
savename = '_setb_2000s_lifetime summary';
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
            if length(dataSave)>1
                if length(dataSave)>2
                    t = 3; 
                else
                    t = 2; 
                end
            else
                t = 1;
            end
            for j  = t %just the second measurement in each set
                datanow = dataSave{j}; 
                curves(count)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2,'color',cm(k,:)); 
                hold all; 
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
savename = '_2000s_degradation';
max_time = 2000; 
time_shift_E = 801610; %amount of time to shift company E measurements over for comparison after switch
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_setb.xlsx'; 
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
control = {'C-1','C-2','FZ','66-2','FZ-12','FZ-new';...
    'Unfired Cz (no H)','Fired Cz (no H)','FZ passivation','mc-Si 750C fired undegraded','FZ degradation','FZ passivation #2'};
fired = {'50b','54b','44b','49b','61b','56b';...
    '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
unfired = {'52b','45b','55b','60b','53b','56b';...
    '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
compE = {'66-2';...
    'mc-Si 750C firing undegraded'};

lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
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
        %Modify the norm calculation for this sample 
        norm_now(:,2) = raw_now(:,2)./raw_now(1,2); 
    end
%     raw_now(:,1) = raw_now(:,1)-time_shift_E; 
%     norm_now(:,1) = norm_now(:,1)-time_shift_E;
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

%% A few helpful plots for analysis
% First plot: initial lifetimes versus MIRHP time
% Second plot: Histogram of measured control lifetimes for mc-Si sample,
% with mean and standard deviation
% Third plot: Histogram of measured control lifetimes for FZ passivation 
% sample with mean and standard deviation
% Fourth plot: Histogram of measured control lifetimes for FZ degradation
% sample with mean and standard deviation
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation';
unfired_names = {'52b','45b','55b','60b'};
unfired_H = [0 10 30 120]; 
fired_names = {'50b','54b','44b','49b'};
fired_H = [0 10 30 120]; 
mcSi_controls = {'53b','61b','56b'}; 
mcSi_controls_H = [30 30 0]; 
mcSi_controls_color = {'r','b','k'}; 
Cz_controls = {'C-1','C-2'}; 
Cz_controls_H = [0 0]; 
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
hgsave(start_tau,[savedirname '\setb_initial_lifetimes']);
print(start_tau,'-dpng','-r0',[savedirname '\setb_initial_lifetimes.png']);

%Second plot
mcSi_tau_controls = {'66-2'}; 
end_times = [0, max_time]; 
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
legend(h,'66-2'); 
set(0,'defaultAxesFontSize', 20)
hgsave(mcSi_controls,[savedirname '\setb_mcSi_histogram']);
print(mcSi_controls,'-dpng','-r0',[savedirname '\setb_mcSi_histogram.png']);

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
hgsave(FZ_pass_control_fig,[savedirname '\setb_FZ_histogram']);
print(FZ_pass_control_fig,'-dpng','-r0',[savedirname '\setb_FZ_histogram.png']);

% %Fourth plot
% FZ_deg_control = {'FZ-12'}; 
% FZ_deg_control_fig = figure('units','normalized','outerposition',[0 0 1 1]);
% for i = 1:length(FZ_deg_control)
%     index = find(strcmp(FZ_deg_control{1,i},samples)==1);
%     norm_now = norm_lifetime_all{index};
%     h(i)=histogram(norm_now(:,2),10); 
%     disp(['Mean for sample ' FZ_deg_control{i} ' is ' num2str(mean(norm_now(:,2)))]); 
%     disp(['Standard deviation for sample ' FZ_deg_control{i} ' is ' num2str(std(norm_now(:,2)))]); 
%     hold on; 
% end
% figure(FZ_deg_control_fig); 
% xlabel('normalized lifetime','FontSize',25); 
% ylabel('frequency','FontSize',25); 
% title('estimating influence of degradation setup on measurement','FontSize',25)
% set(0,'defaultAxesFontSize', 20)
% hgsave(FZ_deg_control_fig,[savedirname '\setb_FZ-12_histogram']);
% print(FZ_deg_control_fig,'-dpng','-r0',[savedirname '\setb_FZ-12_histogram.png']);