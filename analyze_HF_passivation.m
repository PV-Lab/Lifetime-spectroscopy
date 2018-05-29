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
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-28-2018_DegTime624438sec'; 
samples = {'1-5','2-5','3-5','4-5','5-5','6-5','7-5','8-5','P-1'};
% samples = {'P-1'};
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
        meas_res{file,1} = xlsread(this_file,'Summary','N2');%'M2');%'Q2');
        calib{file,1} = xlsread(this_file,'Summary','T2');%'S2');%'T2');
        doping{file,1} = xlsread(this_file,'Summary','E2');
    end
    info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    save([dirname '\' samples{index} '\meas_info.mat'],'info');
end

%% Now analyze the data
clear all; close all; clc;
%Process data after HF passivation

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-28-2018_DegTime624438sec'; 
samples = {'1-5','2-5','3-5','4-5','5-5','6-5','7-5','8-5','P-1'};
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
% savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\Summary\188190s';
% dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\January 11 2018';
% dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\January 15 2018';
% dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\January 17 2018';
% dirnames = {dirname1 dirname2 dirname3}; 
% labels = {'initial','87150s','188190s'};
% savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\Summary\186450s';
% dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\January 11 2018';
% dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\January 16 2018';
% dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\January 18 2018';
% dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#7 degradation\February 1 2018';
% dirnames = {dirname1 dirname2 dirname3 dirname4}; 
% labels = {'initial','initial 2','81480s','186450s'};
% savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\Summary\101735s';
% dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 22 2018';
% dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 23 2018';
% dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 24 2018';
% dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 25 2018';
% dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 26 2018';
% dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 29 2018';
% dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 2 data\#3 degradation\January 31 2018';
% dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7}; 
% labels = {'initial','2500s','5000s','10000s','25000s','50000s','101735s'};
savedirname='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\summary\624438s';
dirname1='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-26-2018_DegTime0sec';
dirname2='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-27-2018_DegTime500sec';
dirname3='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-28-2018_DegTime1000sec';
dirname4='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-29-2018_DegTime2500sec';
dirname5='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-30-2018_DegTime5000sec';
dirname6='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-31-2018_DegTime7500sec';
dirname7='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-02-2018_DegTime10000sec';
dirname8='C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-03-2018_DegTime20000sec';
dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-05-2018_DegTime40000sec';
dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-06-2018_DegTime60000sec';
dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-10-2018_DegTime80000sec'; 
dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-12-2018_DegTime100000sec';
dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-18-2018_DegTime120000sec';
dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-19-2018_DegTime140000sec';
dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-20-2018_DegTime160000sec';
dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-23-2018_DegTime180000sec';
dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-24-2018_DegTime201300sec';
dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-25-2018_DegTime211240sec';
dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-26-2018_DegTime231240sec';
dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\04-30-2018_DegTime248080sec';
dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-02-2018_DegTime266510sec';
dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-03-2018_DegTime284880sec';
dirname23 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-08-2018_DegTime304880sec';
dirname24 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-09-2018_DegTime324880sec';
dirname25 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-10-2018_DegTime348810sec';
dirname26 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-11-2018_DegTime368810sec';
dirname27 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-14-2018_DegTime388810sec';
dirname28 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-15-2018_DegTime408810sec';
dirname29 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-16-2018_DegTime440350sec';
dirname30 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-18-2018_DegTime468875sec';
dirname31 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-21-2018_DegTime500855sec';
dirname32 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-23-2018_DegTime531545sec';
dirname33 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-24-2018_DegTime559338sec';
dirname34 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-25-2018_DegTime598218sec';
dirname35 = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\05-28-2018_DegTime624438sec';
dirnames = {dirname1,dirname2,dirname3,dirname4,dirname5,dirname6,...
    dirname7,dirname8,dirname9,dirname10,dirname11,dirname12,dirname13,...
    dirname14,dirname15,dirname16,dirname17,dirname18,dirname19,dirname20,...
    dirname21,dirname22,dirname23,dirname24,dirname25,dirname26,dirname27,...
    dirname29,dirname30,dirname31,dirname32,dirname33,dirname34,dirname35};
labels= {'0s','500s','1000s','2500s','5000s','7500s','10000s','20000s',...
    '40000s','60000s','80000s','100000s','120000s','140000s','160000s',...
    '180000s','201300s','211240s','231240s','248080s','266510s','284880s',...
    '304880s','324880s','348810s','368810s','388810s','408810s','440350s',...
    '468875','500855','531545s','559338s','598218s','624438s'};
cm = colormap(hsv(length(dirnames))); 
% samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6','P-1'};
% samples = {'1-7','4-7','6-7','7-7'};
% samples = {'2-7','3-7','5-7','8-7'};
% samples = {'1-3','2-3','3-3','4-3','5-3','6-3','7-3','8-3','P-2'};
samples = {'1-5','2-5','3-5','4-5','5-5','6-5','7-5','8-5','P-1'};
savename = '_624438s_lifetime summary';
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\summary\100000s';
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\03-26-2018_DegTime0sec';
% samples = {'1-6','2-6','3-6','4-6','5-6','6-6','7-6','8-6'};
% LTA = [0 600 550 500 450 700 750 650];
samples = {'1-5','2-5','3-5','4-5','5-5','6-5','7-5','8-5'};
LTA = [650 700 450 550 600 0 500 750];
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
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\summary\624438s';
savename = '_624438s_degradation';
max_time = 624438; 
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\PDG\round 3 data\measurement_details.xlsx'; 
deltan_target = 6e14; %target injection level for the measurements, changed to 6e14 on 2/13/17 from 5e14
%Get the measurement details
[meas,samples] = xlsread(meas_details,'measurements');
samples(1,:) = []; 
[times,filenames] = xlsread(meas_details,'filenames'); 
%Make these the same size
filenames = filenames(2:end,2); 

lifetime_all = cell(length(samples),1); 
norm_lifetime_all = cell(length(samples),1); 
cm = colormap(hsv(length(samples))); 
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
% control = {'1-6','5-6','4-6','3-6','2-6','8-6','6-6','7-6','P-1';...
%     'STD','450C','500C','550C','600C','650C','700C','750C','FZ'};
% control = {'1-7','2-7','3-7','4-7','5-7','6-7','7-7','8-7';...
%     'STD','STD','STD','STD','STD','STD','STD','STD',};
% control = {'5-3','7-3','8-3','6-3','4-3','2-3','1-3','3-3','P-2';...
%     'STD','450C','500C','550C','600C','650C','700C','750C','FZ'};
control = {'6-5','3-5','7-5','4-5','5-5','1-5','2-5','8-5','P-1';...
    'STD','450C','500C','550C','600C','650C','700C','750C','FZ'};

LTA = [0;450;500;550;600;650;700;750];
LTA_summary = zeros(length(LTA),5); 
LTA_summary(:,1) = LTA; 

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
    if norm_now(1,1)==0
        norm_now(1,1) = 1; 
    end
    figure(lifetime_raw); 
    loglog(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10,'color',cm(i,:)); 
    hold all; 
    figure(lifetime_norm); 
    semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10,'color',cm(i,:));
    hold all;  
    labels{i,1} = control{2,i};
    %Calculate Nt_star
    Nt_star_now = zeros(length(raw_now),1); 
    for j = 1:length(raw_now)
        Nt_star_now(j) = (1/raw_now(j,2))-(1/raw_now(1,2));  
    end
    figure(Nt_star); 
    loglog(raw_now(:,1),Nt_star_now,'-o','LineWidth',3,'MarkerSize',10,'color',cm(i,:)); 
    hold all; 
    
    max_Nt = find(Nt_star_now == max(Nt_star_now));
    LTA_summary(i,2) = raw_now(max_Nt,1); 
    LTA_summary(i,3) = raw_now(max_Nt,2); 
    LTA_summary(i,4) = norm_now(max_Nt,2);
    LTA_summary(i,5) = Nt_star_now(max_Nt); 
end
figure(lifetime_raw); 
xlabel('time [s]','FontSize',25); 
ylabel('lifetime [s]','FontSize',25); 
legend(labels); 
% title('high hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_raw,[savedirname '\raw' savename]);
print(lifetime_raw,'-dpng','-r0',[savedirname '\raw' savename '.png']);
figure(lifetime_norm); 
xlabel('time [s]','FontSize',25); 
ylabel('norm. lifetime [-]','FontSize',25); 
legend(labels); 
axis([1 max_time 0 2]);
% title('high hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(lifetime_norm,[savedirname '\norm' savename]);
print(lifetime_norm,'-dpng','-r0',[savedirname '\norm' savename '.png']);
figure(Nt_star); 
xlabel('time [s]','FontSize',25); 
ylabel('N_t*','FontSize',25); 
legend(labels,'Location','northwest'); 
% title('high hydrogen','FontSize',25); 
set(0,'defaultAxesFontSize', 20)
hgsave(Nt_star,[savedirname '\Ntstar' savename]);
print(Nt_star,'-dpng','-r0',[savedirname '\Ntstar' savename '.png']);
