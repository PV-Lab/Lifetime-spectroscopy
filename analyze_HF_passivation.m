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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 8 2017\520000s'; 
samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
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

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 8 2017\520000s'; 
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\552550s';
% dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 16 2017';
% dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 17 2017\250s';
% dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 17 2017\500s';
% dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 29 2017\500s';
% dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 29 2017\750s';
% dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 30 2017\1000s';
% dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 31 2017\2000s';
% dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 1 2017\3000s';
% dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 1 2017\4000s';
% dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 5 2017\4000s';
% dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 5 2017\5000s';
% dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 6 2017\6000s';
% dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 7 2017\7000s';
% dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 11 2017\8000s'; 
% dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 14 2017\8000s'; 
% dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 14 2017\9000s'; 
% dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 22 2017\9000s'; 
% dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 22 2017\10000s'; 
% dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 26 2017'; 
% dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 28 2017'; 
% dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 3 2017\20000s'; 
% dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 3 2017\30000s'; 
% dirname23 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 4 2017\40000s';
% dirname24 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 11 2017\40000s'; 
% dirname25 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 12 2017\47000s'; 
% dirname26 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 13 2017\60000s';
% dirname27 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 16 2017\70000s'; 
% dirname28 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 16 2017\80000s';
% dirname29 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 17 2017\90000s';
% dirname30 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 18 2017\100000s'; 
% dirname31 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 20 2017\120000s'; 
% dirname32 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 23 2017\135000s'; 
% dirname33 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 25 2017\150000s';
% dirname34 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 27 2017\155000s'; 
% dirname35 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 30 2017\170000s'; 
% dirname36 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 2 2017\200000s';
% dirname37 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 6 2017\220000s';
% dirname38 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 8 2017\240000s';
% dirname39 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 9 2017\260000s';
% dirname40 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 10 2017\280000s';
% dirnames = {dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 dirname8 dirname9 ...
%     dirname10 dirname11 dirname12 dirname13 dirname14 dirname15 ...
%     dirname16 dirname17 dirname18 dirname19 dirname20 dirname21 ...
%     dirname22 dirname23 dirname24 dirname25 dirname26 dirname27 ...
%     dirname28 dirname29 dirname30 dirname31 dirname32 dirname33 ...
%     dirname34 dirname35 dirname36 dirname37 dirname38 dirname39 ...
%     dirname40}; 
% labels = {'initial','250s','500s','501s','750s','1000s','2000s',...
%     '3000s','4000s','4001s','5000s','6000s','7000s','8000s','9000s',...
%     '9001s','10000s','17300s','20000s','20001s','30000s','40030s','40031s',...
%     '47471s','60221s','70031s','80031s','90031s','100361s','127061s',...
%     '157871s','178661s','190511s','209681s','228941s','255581s',...
%     '285011s','307541s','331001s'};
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 16 2017';
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\August 30 2017\1000s';
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 22 2017\10000s';
dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\September 28 2017';
dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 3 2017\30000s';
dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 4 2017\40000s';
dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 12 2017\47000s';
dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 13 2017\60000s';
dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 16 2017\70000s';
dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 16 2017\80000s';
dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 17 2017\90000s';
dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 18 2017\100000s';
dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 23 2017\135000s';
dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\October 30 2017\170000s';
dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 6 2017\220000s';
dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 9 2017\260000s';
dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 10 2017\280000s';
dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 16 2017\360000s';
dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 17 2017\420000s just after deg';
dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 20 2017\420000s check zeroes';
dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\November 21 2017\450000s';
dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 4 2017\450000s';
dirname23 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 6 2017\490000s';
dirname24 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 8 2017\520000s';
dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
    dirname8 dirname9 dirname10 dirname11 dirname12 dirname13 dirname14 ...
    dirname15 dirname16 dirname18 dirname19 dirname20 dirname21 ...
    dirname22 dirname23 dirname24}; 
labels = {'initial','1000s','10000s','20000s','30000s','40030s',...
    '47471s','60221s','70031s','80031s','90031s','100361s',...
    '157871s','209681s','255581s',...
    '307541s','363310s','421300s just after','421300s wait','456940s',...
    '456940s wait','510190s','552550s'};
cm = colormap(hsv(length(dirnames))); 
samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
savename = '_552550s_lifetime summary';
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
            if length(dataSave)==3
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\552550s';
savename = '_552550s_degradation';
max_time = 552550; 
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
            if length(dataSave)==3
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
    if raw_now(1,1)==0
        raw_now(1,1) = 1; 
    end
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_raw); 
    loglog(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
    hold all;  
    labels{i,1} = control{2,i};
    %Calculate Nt_star
    Nt_star_now = zeros(length(raw_now),1); 
    for j = 1:length(raw_now)
        Nt_star_now(j) = (1/raw_now(j,2))-(1/raw_now(1,2));  
    end
    figure(Nt_star); 
    loglog(raw_now(:,1),Nt_star_now,'-o','LineWidth',3,'MarkerSize',10); 
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
legend(labels,'Location','northwest'); 
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
    if raw_now(1,1)==0
        raw_now(1,1) = 1; 
    end
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_raw); 
    loglog(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
    hold all;  
    labels{i,1} = control{2,i};
    %Calculate Nt_star
    Nt_star_now = zeros(length(raw_now),1); 
    for j = 1:length(raw_now)
        Nt_star_now(j) = (1/raw_now(j,2))-(1/raw_now(1,2));  
    end
    figure(Nt_star); 
    loglog(raw_now(:,1),Nt_star_now,'-o','LineWidth',3,'MarkerSize',10); 
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
legend(labels,'Location','northwest'); 
title('low hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(Nt_star,[savedirname '\lowH_Ntstar' savename]);
print(Nt_star,'-dpng','-r0',[savedirname '\lowH_Ntstar' savename '.png']);

%Make plots to compare each element low vs. high
to_plots = {{'Ti-L-5','Ti-h-5';'Ti low H','Ti high H'},...
    {'Ni-L-5','Ni-h-5';'Ni low H','Ni high H'},...
    {'Mo-L-5','Mo-h-5';'Mo low H','Mo high H'},...
    {'V-L-5','V-h-5';'V low H','V high H'},...
    {'C-L-5','C-h-5';'C low H','C high H'}};
savepre = {'Ti','Ni','Mo','V','C'};

for k = 1:length(to_plots)
    control = to_plots{k}; 
    lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
    Nt_star=figure('units','normalized','outerposition',[0 0 1 1]);
    [nothing,samp] = size(control); 
    labels = {}; 
    for i = 1:samp
        index = find(strcmp(control{1,i},samples)==1);
        raw_now = lifetime_all{index}; 
        if raw_now(1,1)==0
            raw_now(1,1) = 1; 
        end
        norm_now = norm_lifetime_all{index}; 
        figure(lifetime_raw); 
        loglog(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
        hold all; 
        figure(lifetime_norm); 
        semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
        hold all;  
        labels{i,1} = control{2,i};
        %Calculate Nt_star
        Nt_star_now = zeros(length(raw_now),1); 
        for j = 1:length(raw_now)
            Nt_star_now(j) = (1/raw_now(j,2))-(1/raw_now(1,2));  
        end
        figure(Nt_star); 
        loglog(raw_now(:,1),Nt_star_now,'-o','LineWidth',3,'MarkerSize',10); 
        hold all; 
    end
    figure(lifetime_raw); 
    xlabel('time [s]','FontSize',25); 
    ylabel('lifetime [s]','FontSize',25); 
    legend(labels); 
    title('low v. high hydrogen','FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_raw,[savedirname '\lowvhigh_' savepre{k} savename]);
    print(lifetime_raw,'-dpng','-r0',[savedirname '\lowvhigh_' savepre{k} savename '.png']);
    figure(lifetime_norm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('norm. lifetime [-]','FontSize',25); 
    legend(labels); 
    axis([0 max_time 0 2]);
    title('low hydrogen','FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_norm,[savedirname '\lowvhigh_norm_' savepre{k} savename]);
    print(lifetime_norm,'-dpng','-r0',[savedirname '\lowvhigh_norm_' savepre{k} savename '.png']);
    figure(Nt_star); 
    xlabel('time [s]','FontSize',25); 
    ylabel('N_t*','FontSize',25); 
    legend(labels,'Location','northwest'); 
    title('low hydrogen','FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(Nt_star,[savedirname '\lowvhigh_Ntstar_' savepre{k} savename]);
    print(Nt_star,'-dpng','-r0',[savedirname '\lowvhigh_Ntstar_' savepre{k} savename '.png']);
end