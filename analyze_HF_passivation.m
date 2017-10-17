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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 26 2017'; 
samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
% samples = {'Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
for index = 1:length(samples)
%     [fileList,fileListShort] = getAllFiles([dirname '\' samples{index}]); 
    %Get all of the files in this directory
    fileList = dir(fullfile([dirname '\' samples{index}],'*.xlsm')); 
    fileList = {fileList.name}'; 
    savename = [dirname '\' samples{index} '\Raw_data.mat']';
    process_xls_data([dirname '\' samples{index}],savename);
    %We also want to store all of the information for each file
    %T, thickness, resistivity (entered/measured), type, optical constant, calibration,
    %1/64 or 1/1
    for file = 1:length(fileList)
        this_file = fullfile([dirname '\' samples{index}],fileList{file});
%         this_file = fileList{file};
        thick{file,1} = xlsread(this_file,'User','B6');
        res{file,1} = xlsread(this_file,'User','C6');
        oc{file,1} = xlsread(this_file,'User','E6');
        temp{file,1} = 25;
        meas_res{file,1} = xlsread(this_file,'Summary','N2');%'Q2');
        calib{file,1} = xlsread(this_file,'Summary','T2');%'T2');
        doping{file,1} = xlsread(this_file,'Summary','E2');
    end
%     info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    info = struct('filename',fileList,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    save([dirname '\' samples{index} '\meas_info.mat'],'info');
%     clear filelist thick res meas_res oc calib temp doping fileListShort
    clear filelist thick res meas_res oc calib temp doping
end

%% Now analyze the data
clear all; close all; clc;
%Process data after HF passivation

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 26 2017'; 
samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\7000s';
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 16 2017';
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 17 2017\250s';
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 17 2017\500s';
dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 29 2017\500s';
dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 29 2017\750s';
dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 30 2017\1000s';
dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 31 2017\2000s';
dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 1 2017\3000s';
dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 1 2017\4000s';
dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 5 2017\4000s';
dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 5 2017\5000s';
dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 6 2017\6000s';
dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 7 2017\7000s';
dirnames = {dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 dirname8 ...
    dirname10 dirname11 dirname12 dirname13}; 
labels = {'initial','250s','500s','501s','750s','1000s','2000s',...
    '3000s','4000s','4001s','5000s','6000s','7000s'};
cm = colormap(hsv(length(dirnames))); 
samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
savename = '_7000s_lifetime summary';
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\7000s';
savename = '_7000s_degradation';
max_time = 7000; 
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\measurement_details.xlsx'; 
deltan_target = 8e14; %target injection level for the measurements, changed to 6e14 on 2/13/17 from 5e14
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
control = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','17-7-27-1';...
    'Ti','Ni','Mo','V','Control','FZ'};

lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
Nt_star=figure('units','normalized','outerposition',[0 0 1 1]);
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
    %Calculate Nt_star
    Nt_star_now = zeros(length(raw_now),1); 
    for j = 1:length(raw_now)
        Nt_star_now(j) = (1/raw_now(j,2))-(1/raw_now(1,2));  
    end
    figure(Nt_star); 
    plot(raw_now(:,1),Nt_star_now,'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
title('high hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\highH' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\highH' savename '.png']);
figure(lifetime_norm); 
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
axis([0 max_time 0 2]);
title('high hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\highH_norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\highH_norm' savename '.png']);
figure(Nt_star); 
xlabel('time [s]','FontSize',25); 
ylabel('N_t*','FontSize',25); 
legend(labels); 
title('high hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(Nt_star,[savedirname '\highH_Ntstar' savename]);
print(Nt_star,'-dpng','-r0',[savedirname '\highH_Ntstar' savename '.png']);

%Now, which samples do we want to plot together?
control = {'Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-2';...
    'Ti','Ni','Mo','V','Control','FZ'};

lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
Nt_star=figure('units','normalized','outerposition',[0 0 1 1]);
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
    %Calculate Nt_star
    Nt_star_now = zeros(length(raw_now),1); 
    for j = 1:length(raw_now)
        Nt_star_now(j) = (1/raw_now(j,2))-(1/raw_now(1,2));  
    end
    figure(Nt_star); 
    plot(raw_now(:,1),Nt_star_now,'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
title('low hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\lowH' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\lowH' savename '.png']);
figure(lifetime_norm); 
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
axis([0 max_time 0 2]);
title('low hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\lowH_norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\lowH_norm' savename '.png']);
figure(Nt_star); 
xlabel('time [s]','FontSize',25); 
ylabel('N_t*','FontSize',25); 
legend(labels); 
title('low hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(Nt_star,[savedirname '\lowH_Ntstar' savename]);
print(Nt_star,'-dpng','-r0',[savedirname '\lowH_Ntstar' savename '.png']);