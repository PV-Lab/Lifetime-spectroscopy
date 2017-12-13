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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 13 2017'; 
% samples = {'1-7','2-7','3-7','4-7','5-7','6-7','7-7','8-7','66-2','FZ-new',...
%     '17-7-26-STD-1','17-8-7-600-1','17-8-15-550-1','17-8-17-500-1','17-8-31-600-1',...
%     '17-9-1-700-1','17-9-5-450-2','17-9-7-750-1','17-9-11-STD-1'};
% samples = {'17-8-31-650-1',...
%     '17-9-1-700-1','17-9-5-450-2','17-9-7-750-1','17-9-11-STD-1'};
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
        meas_res{file,1} = xlsread(this_file,'Summary','N2');
        calib{file,1} = xlsread(this_file,'Summary','T2');
        doping{file,1} = xlsread(this_file,'Summary','E2');
    end
    info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    save([dirname '\' samples{index} '\meas_info.mat'],'info');
end

%% Now analyze the data
clear all; close all; clc;
%Process data after HF passivation

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 13 2017'; 
% samples = {'1-7','2-7','3-7','4-7','5-7','6-7','7-7','8-7','66-2','FZ-new',...
%     '17-7-26-STD-1','17-8-7-600-1','17-8-15-550-1','17-8-17-500-1','17-8-31-600-1',...
%     '17-9-1-700-1','17-9-5-450-2','17-9-7-750-1','17-9-11-STD-1'};
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\Summary\2000s'; 
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 6 2017'; 
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 7 2017'; 
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 8 2017\before';
dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\December 13 2017';
dirnames = {dirname1 dirname2 dirname3 dirname4}; 
% samples = {'FZ-new',...
%     '17-7-26-STD-1','17-8-7-600-1','17-8-15-550-1','17-8-17-500-1','17-8-31-650-1',...
%     '17-9-1-700-1','17-9-5-450-2','17-9-7-750-1','17-9-11-STD-1'};
% % samples = {'1-7','2-7','3-7','4-7','5-7','6-7','7-7','8-7','66-2'};
samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6','P-1'}; 
times = {'0s','500s','1000s','2000'};
savename = '_2000s_lifetime summary';

curves = [];
label = {};
count = 1;
for i = 1:length(samples)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
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
                t = 2; 
            else
                t = 1;
            end
            datanow = dataSave{t}; 
            curves(count)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2); 
            hold all; 
            label{count} = times{k};
            count = count+1; 
            xlim([5e13 1e17])
        end
    end


xlabel('excess carrier density [cm^-^3]','FontSize',30); 
ylabel('lifetime [s]','FontSize',30);
title(samples{i},'FontSize',30);
legend(curves',label');
set(0,'defaultAxesFontSize', 20)
hgsave(h,[savedirname '\' samples{i} savename]);
print(h,'-dpng','-r0',[savedirname '\' samples{i} savename '.png']);
end

%% Make the degradation curves
clear all; close all; clc; 
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\Degradation\Summary\2000s'; 
savename = '_2000s_lifetime summary';
deltan_target = 8e14; %target injection level for the measurements, changed to 6e14 on 2/13/17 from 5e14
% samples = {'1-7','2-7','3-7','4-7','5-7','6-7','7-7','8-7','66-2','FZ-new',...
%     '17-7-26-STD-1','17-8-7-600-1','17-8-15-550-1','17-8-17-500-1','17-8-31-650-1',...
%     '17-9-1-700-1','17-9-5-450-2','17-9-7-750-1','17-9-11-STD-1'};
samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6','P-1'}; 
lifetime_all = zeros(size(samples)); 

for i = 1:length(samples)
        try 
            load([savedirname '\' samples{i} '\Raw_data.mat']); 
            load([savedirname '\' samples{i} '\meas_info.mat']); 
            flag = 1; 
        catch
            %If that's not successful, just tell me that the sample doesn't
            %exist for that dirname
            warning(['Error accessing data for directory ' num2str(k) ', sample ' samples{i}]);
            flag = 0; 
        end
        if length(dataSave)>1
            t = 2; 
        else
            t = 1;
        end
        datanow = dataSave{t}; 
        [deltan,tau] = remove_duplicates(datanow(:,1),datanow(:,2));
        %There's a little bug in this program, for now just do this... 
        try
            lifetime_all(i) = interp1(deltan,tau,deltan_target); 
        catch
            [deltan,tau] = remove_duplicates(deltan,tau);
            try
                lifetime_all(i) = interp1(deltan,tau,deltan_target);
            catch
                [deltan,tau] = remove_duplicates(deltan,tau);
                lifetime_all(i) = interp1(deltan,tau,deltan_target);
            end
        end
        if isnan(lifetime_all(i))==1
            disp(['Lifetime at ' num2str(deltan_target) ' was NaN for sample ' samples{i}]);
        end
end


%Now, which samples do we want to plot together?
mcsi = {'1-7','2-7','3-7','4-7','5-7','6-7','7-7','8-7','66-2'};
monosi = {'FZ-new',...
    '17-7-26-STD-1','17-8-7-600-1','17-8-15-550-1','17-8-17-500-1','17-8-31-650-1',...
    '17-9-1-700-1','17-9-5-450-2','17-9-7-750-1','17-9-11-STD-1'};

% lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
% [nothing,samp] = size(mcsi); 
% labels = {}; 
% for i = 1:samp
%     index = find(strcmp(control{1,i},samples)==1);
%     raw_now = lifetime_all{index}; 
%     norm_now = norm_lifetime_all{index}; 
%     figure(lifetime_raw); 
%     plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
%     hold all; 
%     figure(lifetime_norm); 
%     plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
%     hold all;  
%     %Get the lifetime, accounting for surface
%     doping_now = control_param(1,i); 
%     W_now = control_param(2,i); 
%     [D_now,Dh] = diffusivity(300,'p',doping_now,deltan_target);
%     tau_rev = zeros(length(raw_now),1); 
%     [num_meas,columns] = size(raw_now); 
%     for j = 1:num_meas
%         try
%             t_index = find(SRV_t==raw_now(j,1)); 
%             tau = raw_now(j,2); 
%             tau_surf = (W_now./(2.*SRV(t_index)))+((1/D_now).*((W_now/pi)^2)); %cm/s
%             tau_rev(j) = ((1./tau)-(1./tau_surf))^(-1); 
%         catch
%             disp(['There was an error calculating the surface lifetime for ' control{1,i} ', time ' num2str(raw_now(j,1)) 's']);
%         end
%     end
%     figure(lifetime_norm_corr); 
%     plot(norm_now(:,1),tau_rev./tau_rev(1),'-o','LineWidth',3,'MarkerSize',10);
%     hold all; 
%     labels{i,1} = control{2,i};
% end
% figure(lifetime_raw); 
% xlabel('time [s]','FontSize',25); 
% ylabel('lifetime [s]','FontSize',25); 
% legend(labels); 
% title('control samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_raw,[savedirname '\controls' savename]);
% print(lifetime_raw,'-dpng','-r0',[savedirname '\controls' savename '.png']);
% figure(lifetime_norm); 
% xlabel('time [s]','FontSize',25); 
% ylabel('norm. lifetime [-]','FontSize',25); 
% legend(labels); 
% axis([0 max_time 0 2]);
% title('control samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_norm,[savedirname '\controls_norm' savename]);
% print(lifetime_norm,'-dpng','-r0',[savedirname '\controls_norm' savename '.png']);
% figure(lifetime_norm_corr); 
% xlabel('time [s]','FontSize',25); 
% ylabel('norm. lifetime (no surface) [-]','FontSize',25); 
% legend(labels); 
% axis([0 max_time 0 2]);
% title('control samples (surface removed)','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_norm_corr,[savedirname '\controls_norm_nosurf' savename]);
% print(lifetime_norm_corr,'-dpng','-r0',[savedirname '\controls_norm_nosurf' savename '.png']);
%     
% lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
% lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
% [nothing,samp] = size(unfired); 
% labels = {};
% for i = 1:samp
%     index = find(strcmp(unfired{1,i},samples)==1);
%     raw_now = lifetime_all{index}; 
%     norm_now = norm_lifetime_all{index}; 
%     figure(lifetime_raw); 
%     plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
%     hold all; 
%     figure(lifetime_norm); 
%     plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
%     hold all; 
%     labels{i,1} = unfired{2,i}; 
% end
% figure(lifetime_raw); 
% xlabel('time [s]','FontSize',25); 
% ylabel('lifetime [s]','FontSize',25); 
% legend(labels); 
% title('unfired samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_raw,[savedirname '\unfired' savename]);
% print(lifetime_raw,'-dpng','-r0',[savedirname '\unfired' savename '.png']);
% figure(lifetime_norm); 
% xlabel('time [s]','FontSize',25); 
% ylabel('norm. lifetime [-]','FontSize',25); 
% legend(labels); 
% axis([0 max_time 0 2]);
% title('unfired samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_norm,[savedirname '\unfired_norm' savename]);
% print(lifetime_norm,'-dpng','-r0',[savedirname '\unfired_norm' savename '.png']);
% 
% lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
% lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
% [nothing,samp] = size(fired); 
% labels = {};
% for i = 1:samp
%     index = find(strcmp(fired{1,i},samples)==1);
%     raw_now = lifetime_all{index}; 
%     norm_now = norm_lifetime_all{index}; 
%     figure(lifetime_raw); 
%     plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
%     hold all; 
%     figure(lifetime_norm); 
%     plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
%     hold all; 
%     labels{i,1} = fired{2,i}; 
% end
% figure(lifetime_raw); 
% xlabel('time [s]','FontSize',25); 
% ylabel('lifetime [s]','FontSize',25); 
% legend(labels); 
% title('fired samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_raw,[savedirname '\fired' savename]);
% print(lifetime_raw,'-dpng','-r0',[savedirname '\fired' savename '.png']);
% figure(lifetime_norm); 
% axis([0 max_time 0 2]);
% xlabel('time [s]','FontSize',25); 
% ylabel('norm. lifetime [-]','FontSize',25); 
% legend(labels); 
% title('fired samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_norm,[savedirname '\fired_norm' savename]);
% print(lifetime_norm,'-dpng','-r0',[savedirname '\fired_norm' savename '.png']);
% 
% lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
% lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
% [nothing,samp] = size(compE); 
% labels = {};
% for i = 1:samp
%     index = find(strcmp(compE{1,i},samples)==1);
%     raw_now = lifetime_all{index}; 
%     norm_now = norm_lifetime_all{index};
%     if strcmp(compE{1,i},'68-2')==1 
%         remove = find(raw_now(:,1)<time_shift_E); 
%         raw_now(remove,:) = []; 
%         remove = find(norm_now(:,1)<time_shift_E); 
%         norm_now(remove,:) = []; 
%         %Modify the norm calculation for this sample 
%         norm_now(:,2) = raw_now(:,2)./raw_now(1,2); 
%     end
%     raw_now(:,1) = raw_now(:,1)-time_shift_E; 
%     norm_now(:,1) = norm_now(:,1)-time_shift_E;
%     figure(lifetime_raw); 
%     plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
%     hold all; 
%     figure(lifetime_norm); 
%     plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
%     hold all; 
%     labels{i,1} = compE{2,i}; 
% end
% figure(lifetime_raw); 
% xlabel('time [s]','FontSize',25); 
% ylabel('lifetime [s]','FontSize',25); 
% legend(labels); 
% title('Company E samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_raw,[savedirname '\compE' savename]);
% print(lifetime_raw,'-dpng','-r0',[savedirname '\compE' savename '.png']);
% figure(lifetime_norm); 
% axis([0 max_time 0 2]);
% xlabel('time [s]','FontSize',25); 
% ylabel('norm. lifetime [-]','FontSize',25); 
% legend(labels); 
% title('Company E samples','FontSize',25); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(lifetime_norm,[savedirname '\compE_norm' savename]);
% print(lifetime_norm,'-dpng','-r0',[savedirname '\compE_norm' savename '.png']);

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