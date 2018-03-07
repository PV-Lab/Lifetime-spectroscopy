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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 9 2018\1350000s'; 
% samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
samples = {'C-L-5'};
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

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 9 2018\1350000s'; 
% samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
samples = {'C-L-5'};
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
%     try
%         lifetime_store(i) = interp1(deltan,tau,1e15); 
%     catch
%         [deltan,tau] = remove_duplicates(deltan,tau);
%         try 
%             lifetime_store(i) = interp1(deltan,tau,1e15);
%         catch
%             [deltan,tau] = remove_duplicates(deltan,tau);
%             lifetime_store(i) = interp1(deltan,tau,1e15);
%         end
%     end
end

%% Analyze different states together
clear all; close all; clc;
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\1814530s';
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
dirname25 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 15 2017\540000s';
dirname26 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 18 2017\570000s';
dirname27 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\December 18 2017\just after deg';
dirname28 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 5 2018\595000s'; 
dirname29 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 8 2018\620000s';
dirname30 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 9 2018\670000s';
dirname31 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 11 2018\740000s';
dirname32 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 15 2018\820000s'; 
dirname33 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 17 2018\900000s';
dirname34 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 22 2018\1000000s';
dirname35 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\January 31 2018\1100000s';
dirname36 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 5 2018\1250000s';
dirname37 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 9 2018\1350000s';
dirname38 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 14 2018\1350000s';
dirname39 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 19 2018\1530000s';
dirname40 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\February 23 2018\1650000s';
dirname41 = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\March 6 2018';
dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
    dirname8 dirname9 dirname10 dirname11 dirname12 dirname13 dirname14 ...
    dirname15 dirname16 dirname18 dirname19 dirname20 dirname21 ...
    dirname22 dirname23 dirname24 dirname25 dirname26 dirname27 ...
    dirname28 dirname29 dirname30 dirname31 dirname32 dirname33 ...
    dirname34 dirname35 dirname36 dirname37 dirname38 dirname39 ...
    dirname40 dirname41}; 
labels = {'initial','1000s','10000s','20000s','30000s','40030s',...
    '47471s','60221s','70031s','80031s','90031s','100361s',...
    '157871s','209681s','255581s',...
    '307541s','363310s','421300s just after','421300s wait','456940s',...
    '456940s wait','510190s','552550s','599290s','648910s','674200s just after',...
    '674200s wait','699070s','757210s','842440s','937360s','1024120s',...
    '1106680s','1202380s','1357600s','1468120s','1468120s wait','1648300s',...
    '1814530s','1814530s wait'};
cm = colormap(hsv(length(dirnames))); 
samples = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-1','17-7-27-2'};
savename = '_1814530s_lifetime summary';
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\1814530s';
savename = '_1814530s_degradation';
max_time = 1814530; 
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

%% Analyze lifetime throughout degradation
%Things we want to do: 
%1) Load injection dependent data for each sample
%2) determine the SRV from the FZ samples (average?)
%3) Apply SRV to each sample to get SRH lifetime
%4) Fit SRH lifetime with 1 or 2 defects? Select the dominant defect
%5) Plot the k value from midgap for each sample throughout degradation 
%Key comparisons: directly after degradation, versus after a wait
%Issues remaining: Do we need to stitch for
%the "C" samples? Need to figure out which is the dominant defect. --
%should we choose just a limited number of measurements to analyze?
clear all; close all; clc; 
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\1814530s\lifetime spectroscopy';
savename = '_1814530s_degradation';
max_time = 1814530; 
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\measurement_details.xlsx'; 
deltan_target = 8e14; %target injection level for the measurements, changed to 6e14 on 2/13/17 from 5e14
%Get the measurement details
[meas,samples] = xlsread(meas_details,'measurements');
samples(1,:) = []; 
[times,filenames] = xlsread(meas_details,'filenames'); 
%Make these the same size
filenames = filenames(2:end,2); 
temp = 25; %degrees C
type = 'p'; 

%Which are the float-zone samples?
floatzone = {'17-7-27-1','17-7-27-2'};

all_data = cell(length(samples),3); 
    
for i = 1:length(samples)
    meas_thissample = meas(i,:);
    data_thissample = cell(length(meas_thissample),3); 
    for j = 1:length(meas_thissample)
        if isnan(meas_thissample(j))==0
            %save the time for this sample
            data_thissample{j,1} = meas_thissample(j); 
            %now create the proper filename
            findex = find(meas_thissample(j)==times);  
            filename = [filenames{findex} '\' samples{i} '\Raw_data.mat'];
            load(filename);
            %the order of files is typically:
            %1 - 1/64 (if measured)
            %2 - avg3 (normal height)
            %3 - high_avg3 (high height)
            if length(dataSave)==3
                t = 3; %choose the high height measurement (C)
            elseif length(dataSave)==2
                t = 1; %choose the normal height measurement (mcSi)
            else  
                t = 1; %this will be for the FZ samples where there's only one measurement
            end
            datanow = dataSave{t}; 
            [deltan,tau] = remove_duplicates(datanow(:,1),datanow(:,2));
            %save the injection-dependent data
            data_thissample{j,2} = deltan; 
            data_thissample{j,3} = tau; 
            if j == 1
                %grab the doping and the thickness
                filename = [filenames{findex} '\' samples{i} '\meas_info.mat'];
                load(filename);
                all_data{i,2} = info(t).thickness;
                all_data{i,3} = info(t).doping; 
            end
        end
    end
    all_data{i,1} = data_thissample; 
end

%Now we evaluate the SRV at each degradation point
SRV_store = cell(size(floatzone)); 
for i = 1:length(floatzone)
    index = find(strcmp(samples,floatzone{i})==1);
    all_data_now = all_data{index,1}; 
    [num_meas,n] = size(all_data_now); 
    SRV_thissample = cell(num_meas,3); 
    for k = 1:num_meas
        deltan = all_data_now{k,2}; 
        tau = all_data_now{k,3}; 
        tau_intr = zeros(size(deltan));
        diffusivity_save = zeros(size(deltan)); 
        for j = 1:length(deltan)
            %Get the intrinsic lifetime
            tau_intr(j,1) = Richter(temp+273.15,deltan(j),all_data{index,3},type);
            %unfortunately the diffusivity is also temperature and injection
            %dependent
            [De,Dh] = diffusivity(temp+273.15,type,all_data{index,3},deltan(j));
            if type == 'n'
                diffusivity_save(j,1) = Dh; %hole is minority carrier
            elseif type == 'p'
                diffusivity_save(j,1) = De; %electron is minority carrier
            end
        end
        %Calculate the surface-related lifetime assuming zero SRH contribution
        tau_surf = ((1./tau)-(1./tau_intr)).^(-1);
        %Calculate the SRV, including the injection-dependent diffusivity
        SRV = all_data{index,2}./((tau_surf-((1./diffusivity_save).*((all_data{index,2}/pi)^2))).*2);
        SRV_thissample{k,3} = SRV; 
        SRV_thissample{k,2} = deltan; 
        SRV_thissample{k,1} = all_data_now{k,1}; 
    end
    SRV_store{i} = SRV_thissample; 
end

%Go through and find an average for the surface contribution to lifetime
baseFZ = SRV_store{1}; 
SRV_avg = cell(size(baseFZ)); 
[num_meas,n] = size(SRV_avg); 
figure;
label = {};
for i = 1:num_meas
    SRV_sum = baseFZ{i,3}; 
    base_deltan = baseFZ{i,2}; 
    count = 1;
    for k = 2:length(floatzone)
        thisFZ = SRV_store{k}; 
        index = find(strcmp(baseFZ{i,1},thisFZ(:,1))==1); 
        if isempty(index)==0
            try
                SRV_now = interp1(thisFZ{index,2},thisFZ{index,3},base_deltan); 
            catch
                [thisFZ{index,2},thisFZ{index,3}] = remove_duplicates(thisFZ{index,2},thisFZ{index,3});
                try
                    SRV_now = interp1(thisFZ{index,2},thisFZ{index,3},base_deltan); 
                catch
                    [thisFZ{index,2},thisFZ{index,3}] = remove_duplicates(thisFZ{index,2},thisFZ{index,3});
                    SRV_now = interp1(thisFZ{index,2},thisFZ{index,3},base_deltan); 
                end
            end
            SRV_sum = SRV_sum+SRV_now; 
            count = count+1; 
        end
    end
    SRV_avg{i,3} = SRV_sum/count; 
    SRV_avg{i,1} = baseFZ{i,1}; 
    SRV_avg{i,2} = baseFZ{i,2}; 
    semilogx(base_deltan,SRV_avg{i,3},'LineWidth',2); 
    label{i,1} = num2str(baseFZ{i,1}); 
    hold all;
end
xlabel('excess carrier density [cm^-^3]','FontSize',20); 
ylabel('SRV [cm/s]','FontSize',20); 
legend(label); 

lifetime_breakdown = figure; 
defect_parameters = cell(length(samples)-length(floatzone),2); 
%Now, for each sample, calculate and then fit the SRH lifetime
for i = 1:length(samples)
    %Check to make sure it's not a FZ sample
    if strcmp(samples{i},floatzone(:))==0 
        %set the low/high injection cutoffs
        if strcmp(samples{i},'C-h-5')==1 || strcmp(samples{i},'C-L-5')==1
            cutoff_low = 5e13;
            cutoff_high = 5e15; 
        else
            cutoff_low = 5e13;
            cutoff_high = 8e15; 
        end
        linSRH = figure; 
        label = {}; 
        all_data_now = all_data{i,1}; 
        [num_meas,n] = size(all_data_now); 
        def_param_thissample = cell(num_meas,5);
        easy_summary = zeros(num_meas,6); 
        for j = 1:num_meas
            SRV_index = find(all_data_now{j,1}==cell2mat(SRV_avg(:,1)));
            SRV = SRV_avg{SRV_index,3};
            deltanSRV = SRV_avg{SRV_index,2}; 
            deltan = all_data_now{j,2}; 
            tau = all_data_now{j,3}; 
            tau_intr = zeros(size(deltan));
            diffusivity_save = zeros(size(deltan)); 
            for k = 1:length(deltan)
                %Get the intrinsic lifetime
                tau_intr(k,1) = Richter(temp+273.15,deltan(k),all_data{i,3},type);
                %unfortunately the diffusivity is also temperature and injection
                %dependent
                [De,Dh] = diffusivity(temp+273.15,type,all_data{i,3},deltan(k));
                if type == 'n'
                    diffusivity_save(k,1) = Dh; %hole is minority carrier
                elseif type == 'p'
                    diffusivity_save(k,1) = De; %electron is minority carrier
                end
            end
            SRVq = interp1(deltanSRV,SRV,deltan); 
            tau_surf =(all_data{i,2}./(2.*SRVq))+((1./diffusivity_save(1)).*((all_data{i,2}/pi)^2)); %cm/s
            %Finally, we can calculate the SRH lifetime
            tau_SRH = ((1./tau)-(1./tau_intr)-(1./tau_surf)).^(-1);
            indices = find(tau_SRH<0 | isnan(tau_SRH)==1); 
            tau_SRH(indices) = [];tau(indices) = []; tau_intr(indices) = [];
            tau_surf(indices) = []; deltan(indices) = [];
            %Plot the different mechanisms and ask to crop each one
            figure(lifetime_breakdown); clf; 
            loglog(deltan,tau); hold all; loglog(deltan,tau_surf); 
            hold all; loglog(deltan,tau_intr); hold all; 
            loglog(deltan,tau_SRH); 
%             xlabel('excess carrier density [cm^-^3]'); 
%             ylabel('lifetime [s]'); 
%             legend('measured','surface','intrinsic','SRH'); 
%             disp('Select the region for cutting off the HIGH injection data');
%             [cutoff,nothing]=ginput(1);
            [deltan,tau_SRH] = remove_highinj(deltan,tau_SRH,cutoff_high);
%             disp('Select the region for cutting off the LOW injection data');
%             [cutoff,nothing]=ginput(1);
            [deltan,tau_SRH] = remove_lowinj(deltan,tau_SRH,cutoff_low);
            hold all;
            loglog(deltan,tau_SRH,'x'); 
            hgsave(lifetime_breakdown,[savedirname '\' samples{i} '_' num2str(all_data_now{j,1}) 's_lifetimeBreakdown']);
            %Now that we have this, linearize and fit
            [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(temp+273.15,all_data{i,3},type); 
            %Normalized carrier density
            if type == 'p'
                X = (n0+deltan)./(p0+deltan);
            elseif type == 'n'
                X = (p0+deltan)./(n0+deltan);
            end
            figure(linSRH);
            plot(X,tau_SRH,'LineWidth',2); 
            hold all; 
            label{j,1} = num2str(all_data_now{j,1}); 
            %We assume (based on some previous analysis) that we have two
            %defects present at each measured point
            savename = [savedirname '\' samples{i} '_' num2str(all_data_now{j,1}) 's_'];
            [two_defects,MSE_two,all_parameters_store,all_MSE_store] = fit_murphy_two(X,tau_SRH,temp,savename,1e5);
            def_param_thissample{j,1} = two_defects; 
            def_param_thissample{j,2} = MSE_two; 
            [m,n] = size(two_defects);
            Et = cell(size(two_defects)); k = cell(size(two_defects)); 
            alphanN = cell(size(two_defects)); 
            for x = 1:m
                [Et{x},k{x},alphanN{x}]=generate_Ek(two_defects(x,:),temp+273.15,all_data{i,3},type);
            end
            def_param_thissample{j,3} = Et; 
            def_param_thissample{j,4} = k; 
            def_param_thissample{j,5} = alphanN; 
            %Pick the dominant defect in low injection (3e14 to be safe)
            inj = find(abs(deltan-8e14)==min(abs(deltan-8e14))); 
            actual = tau_SRH(inj); 
            def1 = two_defects(1,1)*X(inj)+two_defects(2,2);
            def2 = two_defects(2,1)*X(inj)+two_defects(2,2);
%             Eg = Sze(temp+273.15); 
            %Revise the below for defect 2 - save data!! 
            if abs(def1-actual)<abs(def2-actual)
                Et_now = Et{1}; 
                Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
                k_now = k{1}; 
                alphanN_now = alphanN{1}; 
                %Defect 1 is dominant
                to_write = [all_data_now{j,1} k_now(Et_index) alphanN_now(Et_index)];
                %Defect 2 is secondary
                Et_now = Et{2}; 
                Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
                k_now = k{2}; 
                alphanN_now = alphanN{2};
                to_write(1,4) = k_now(Et_index);
                to_write(1,5) = alphanN_now(Et_index);
                to_write(1,6) = MSE_two; 
                easy_summary(j,:) = to_write; 
            elseif abs(def2-actual)<abs(def1-actual)
                Et_now = Et{2}; 
                Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
                k_now = k{2}; 
                alphanN_now = alphanN{2}; 
                %Defect 2 is dominant
                to_write = [all_data_now{j,1} k_now(Et_index) alphanN_now(Et_index)];
                %Defect 1 is secondary
                Et_now = Et{1}; 
                Et_index = find(abs(0-Et_now)==min(abs(0-Et_now))); 
                k_now = k{1}; 
                alphanN_now = alphanN{1};
                to_write(1,4) = k_now(Et_index);
                to_write(1,5) = alphanN_now(Et_index);
                to_write(1,6) = MSE_two; 
                easy_summary(j,:) = to_write;
           end
        end
        defect_parameters{i,1} = def_param_thissample; 
        defect_parameters{i,2} = easy_summary; 
        figure(linSRH);
        xlabel('X [-]','FontSize',20); 
        ylabel('tau_S_R_H [s]','FontSize',20);
        legend(label); 
        title(samples{i},'FontSize',20); 
        hgsave(linSRH,[savedirname '\' samples{i} '_linSRH']);
        print(linSRH,'-dpng','-r0',[savedirname '\' samples{i} '_linSRH.png']);
    end
end

%Save all of the important data
save([savedirname '\lifetime_spect_analysis'],'all_data','SRV_store',...
    'SRV_avg','defect_parameters','type','temp','max_time'); 

%% Plot the k-values and lifetime degradation on comparable figure
%Subplot 3x2. First column = low hydrogen. Second column = high hydrogen. 
%First row = lifetime, or normalized lifetime, or effect defect density
%throughout degradation?
%Second row = k value defect 1 at the same time points throughout degradation
%Third row = k value defect 2 at the same time points throughout
%degradation

clear all; close all; clc; 
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\Summary\1814530s\lifetime spectroscopy';
savename = '_1814530s_degradation';
max_time = 1814530; 
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
            %the order of files is typically:
            %1 - 1/64 (if measured)
            %2 - avg3 (normal height)
            %3 - high_avg3 (high height)
            if length(dataSave)==3
                t = 3; %choose the high height measurement (C)
            elseif length(dataSave) == 2
                t = 1; %choose the normal height measurement (mcSi)
            else  
                t = 1; %this will be for the FZ samples where there's only one measurement
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

%Now that we have the lifetime data - normalized, raw - load the fitted
%lifetime data from the previous section. 
load([savedirname '\lifetime_spect_analysis.mat']); 
    
%Now, which samples do we want to plot together?
control = {'Ti-L-5','Ni-L-5','Mo-L-5','V-L-5','C-L-5','17-7-27-2';...
    'Ti','Ni','Mo','V','Control','FZ'};

%We will have all the samples on the same figure. 
lifetime_norm = figure; 
[nothing,samp] = size(control); 
labels = {}; 
for i = 1:samp
    %First plot the normalized lifetime for this sample
    index = find(strcmp(control{1,i},samples)==1);
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_norm); 
    subplot(4,2,1); 
    if norm_now(1,1) == 0
        norm_now(1,1) = 1; 
    end
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',2,'MarkerSize',5); 
    hold all; 
    labels{i,1} = control{2,i};
    %only do the following if we are not working with FZ
    if strcmp(control{2,i},'FZ')==0
        %Now we get the defect parameters which have been processed in the
        %previous section
        defects_now = defect_parameters{index,2}; %easy summary
        if defects_now(1,1) == 0
            defects_now(1,1) = 1; 
        end
        %defect 1
        subplot(4,2,3); 
        loglog(defects_now(:,1),defects_now(:,2),'-o','LineWidth',2,'MarkerSize',5);
        hold all;
        %defect 2
        subplot(4,2,5); 
        loglog(defects_now(:,1),defects_now(:,4),'-o','LineWidth',2,'MarkerSize',5); 
        hold all; 
        %MSE
        subplot(4,2,7); 
        loglog(defects_now(:,1),defects_now(:,6),'-o','LineWidth',2,'MarkerSize',5); 
        hold all; 
    end
end

%Lets format these plots
subplot(4,2,1); 
% xlabel('degradation time [s]','FontSize',14); 
ylabel('norm. lifetime [-]','FontSize',14); 
title('low hydrogen','FontSize',14); 
legend(labels); 
axis([1 max_time*1.25 0 2]);
subplot(4,2,3); 
% xlabel('degradation time [s]','FontSize',14); 
ylabel('k-value [-]','FontSize',14);  
title('dominant defect','FontSize',14); 
axis([1 max_time*1.25 0.4 100]);
subplot(4,2,5); 
% xlabel('degradation time [s]','FontSize',14); 
ylabel('k-value [-]','FontSize',14);  
title('secondary defect','FontSize',14); 
axis([1 max_time*1.25 1 1000]);
subplot(4,2,7); 
xlabel('degradation time [s]','FontSize',14); 
ylabel('MSE fit','FontSize',14);  
xlim([1 max_time*1.25]);

%Do the same for the other hydrogen samples
%Now, which samples do we want to plot together?
control = {'Ti-h-5','Ni-h-5','Mo-h-5','V-h-5','C-h-5','17-7-27-1';...
    'Ti','Ni','Mo','V','Control','FZ'};
 
[nothing,samp] = size(control); 
labels = {}; 
for i = 1:samp
    %First plot the normalized lifetime for this sample
    index = find(strcmp(control{1,i},samples)==1);
    norm_now = norm_lifetime_all{index}; 
    figure(lifetime_norm); 
    subplot(4,2,2); 
    if norm_now(1,1) == 0
        norm_now(1,1) = 1; 
    end
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',2,'MarkerSize',5); 
    hold all; 
    labels{i,1} = control{2,i};
    %only do the following if we are not working with FZ
    if strcmp(control{2,i},'FZ')==0
        %Now we get the defect parameters which have been processed in the
        %previous section
        defects_now = defect_parameters{index,2}; %easy summary
        if defects_now(1,1) == 0
            defects_now(1,1) = 1; 
        end
        %defect 1
        subplot(4,2,4); 
        loglog(defects_now(:,1),defects_now(:,2),'-o','LineWidth',2,'MarkerSize',5);
        hold all;
        %defect 2
        subplot(4,2,6); 
        loglog(defects_now(:,1),defects_now(:,4),'-o','LineWidth',2,'MarkerSize',5); 
        hold all;
        %MSE
        subplot(4,2,8); 
        loglog(defects_now(:,1),defects_now(:,6),'-o','LineWidth',2,'MarkerSize',5); 
        hold all; 
    end
end

%Lets format these plots
subplot(4,2,2); 
% xlabel('degradation time [s]','FontSize',14); 
ylabel('norm. lifetime [-]','FontSize',14); 
title('high hydrogen','FontSize',14); 
legend(labels); 
axis([1 max_time*1.25 0 2]);
subplot(4,2,4); 
% xlabel('degradation time [s]','FontSize',14); 
ylabel('k-value [-]','FontSize',14);  
title('dominant defect','FontSize',14); 
axis([1 max_time*1.25 0.4 100]);
subplot(4,2,6); 
% xlabel('degradation time [s]','FontSize',14); 
ylabel('k-value [-]','FontSize',14);  
title('secondary defect','FontSize',14); 
axis([1 max_time*1.25 1 1000]);
subplot(4,2,8); 
xlabel('degradation time [s]','FontSize',14); 
ylabel('MSE fit','FontSize',14);  
xlim([1 max_time*1.25]);
hgsave(lifetime_norm,[savedirname '\fitted_lifetime_summary']);
print(lifetime_norm,'-dpng','-r0',[savedirname '\fitted_lifetime_summary.png']);