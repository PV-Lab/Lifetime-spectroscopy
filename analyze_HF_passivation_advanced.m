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
%This script performs analysis for samples degraded and measured with HF
%passivation. This includes improvements for switching between different
%sample sets as well as for analyzing the lifetime independent of surface
%fluctuations. 

%% Run this first to set the ground rules. Edit the function HFbora with new information
clear all; close all; clc; 
%Change these values
%-----------------------------
bora = 'compE'; %'set-b' or 'set-a' or 'compE'
%Most recent directory that we want to analyze now. 
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 24 2017';
%where we want to save any new, non-sample-specific data
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation'; 
%Spreadsheet specification for the actual measurements
spreadsheet = 'new'; %old (before TS) or new (after TS)
%target injection level for the measurements, used to make degradation
%curves
deltan_target = 6e14; %cm^-3

%-----------------------------
%which sample set, b or a? 
[samples,dirnames,labels,savename,surface_control,...
    plotting_group,plotting_names,meas_details,max_time]=HFbora(bora);

%% First process the raw data
for index = 1:length(samples)
    [fileList,fileListShort] = getAllFiles([dirname '\' samples{index}]); 
    savename_excel = [dirname '\' samples{index} '\Raw_data.mat']';
    try
        process_xls_data([dirname '\' samples{index}],savename_excel);
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
            doping{file,1} = xlsread(this_file,'Summary','E2');
            temp{file,1} = 25;
            meas_res{file,1} = xlsread(this_file,'Summary','N2');
            if strcmp(spreadsheet,'new')==1
                calib{file,1} = xlsread(this_file,'Summary','T2');
            elseif strpcmp(spreadsheet,'old')==1
                calib{file,1} = xlsread(this_file,'Summary','S2');
            else
                disp('There was an error defining spreadsheet old or new'); 
            end
        end
        info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
        save([dirname '\' samples{index} '\meas_info.mat'],'info');
    catch 
        disp(['Error reading sample ' samples{index}]); 
    end
end
%% Now analyze the data and plot lifetime curves for this measurement
close all; clc; 
lifetime_store = zeros(length(samples),1); 
for i = 1:length(samples)
    try
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
    catch 
        disp(['Error loading data on this date for sample ' samples{i}]); 
    end
end
%% Now we want to plot lifetime curves together across different deg times
close all; clc; 
cm = colormap(hsv(length(dirnames)));
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
    hgsave(h,[savedirname '\' samples{i} savename '_lifetime summary']);
    print(h,'-dpng','-r0',[savedirname '\' samples{i} savename '_lifetime summary.png']);
end
%% Make the degradation curves
close all; clc; 
%EIGHT types: 
% 1) raw lifetime at a specific injection level throughout degradation
% (lifetime_all)
% 2) normalized lifetime at a specific injection level throughout
% degradation (norm_lifetime_all)
% 3) lifetime at a specific injection level, corrected by the surface
% lifetime from the FZ sample (FZ_corr)
% 4) #3, normalized (FZ_corr_norm)
% 5) lifetime at a specific injection level, corrected by the surface
% lifetime from the mc-Si sample using harmonic difference (mcSi_harm)
% 6) #5, normalized (mcSi_harm_norm)
% 7) lifetime at a specific injection level, corrected by the surface
% lifetime from the mc-Si sample using a ratio to this value (mcSi_ratio)
% 8) #7, normalized (mcSi_ratio_norm)

%Get the measurement details
[meas,samples] = xlsread(meas_details,'measurements');
samples(1,:) = []; 
[times,filenames] = xlsread(meas_details,'filenames'); 
%Make these the same size
filenames = filenames(2:end,2); 
%Initialize the structures
lifetime_all = cell(length(samples),1); 
norm_lifetime_all = cell(length(samples),1); 
thickness_all = cell(length(samples),1); 
doping_all = cell(length(samples),1); 
for i = 1:length(samples)
    meas_thissample = meas(i,:);
    lifetime_store = [];
    doping_store = [];
    thickness_store = []; 
    for j = 1:length(meas_thissample)
        if isnan(meas_thissample(j))==0
            %now create the proper filename
            findex = find(meas_thissample(j)==times);  
            filename = [filenames{findex} '\' samples{i} '\Raw_data.mat'];
            load(filename);
            filename = [filenames{findex} '\' samples{i} '\meas_info.mat'];
            load(filename);
            if length(dataSave)>1
                t = 2; 
            else
                t = 1;
            end
            doping_store(j) = info(t).doping; 
            thickness_store(j) = info(t).thickness; 
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
    thickness_all{i} = thickness_store; 
    doping_all{i} = doping_store; 
end

%Now we need to make the corrections - first, using FZ as the surface
%reference.
index_FZ = [];
for i = 1:length(surface_control)
    index_FZ(end+1) = find(strcmp(samples,surface_control{i})==1); 
end
%We need to get the surface information from FZ before we loop
SRV_FZ = cell(size(index_FZ)); 
SRV_label = cell(size(index_FZ)); 
figure; 
for i = 1:length(index_FZ)
    doping_FZ = doping_all{index_FZ(i)}; %cm-3
    thickness_FZ = thickness_all{index_FZ(i)}; %cm
    %Check that the parameters are consistent for this sample
    if length(doping_FZ) ~= length(find(doping_FZ(1)==doping_FZ))
        disp(['the doping for sample ' samples{index_FZ(i)} ' is not consistent']); 
    end
    if length(thickness_FZ) ~= length(find(thickness_FZ(1)==thickness_FZ))
        disp(['the thickness for sample ' samples{index_FZ(i)} ' is not consistent']); 
    end
    raw_FZ = lifetime_all{index_FZ(i)}; 
    [meas,nothing]=size(raw_FZ); 
    SRV = zeros(meas,1); 
    for j = 1:length(SRV)
        [De,Dh] = diffusivity(300,'p',doping_FZ(j),deltan_target); 
        D_FZ  = De; %cm2/s
        tau_intr = Richter(300,deltan_target,doping_FZ(j),'p');
        tau_surf = ((1./raw_FZ(j,2))-(1./tau_intr)).^(-1);
        SRV(j) = thickness_FZ(j)./((tau_surf-((1/D_FZ)*((thickness_FZ(j)/pi)^2))).*2);
    end
    SRV_FZ{i} = [raw_FZ(:,1) SRV]; 
    plot(raw_FZ(:,1),SRV); 
    hold all; 
    SRV_label{i} = samples{index_FZ(i)}; 
end
axis([0 max_time 0 10]);
xlabel('time [s]','FontSize',25); 
ylabel('SRV [cm/s]','FontSize',25); 
legend(SRV_label); 

FZ_corr = cell(size(samples)); 
FZ_corr_norm = cell(size(samples)); 
count = 1; 
for i = 1:length(samples)
    if strcmp(samples{i},'FZ')==0 && strcmp(samples{i},'FZ-new')==0 && ...
            strcmp(samples{i},'68-4')==0 && strcmp(samples{i},'60a')==0 && ...
            strcmp(samples{i},'56b')==0
        raw_now = lifetime_all{i};
        doping_now = doping_all{i}; 
        thickness_now = thickness_all{i}; 
        %Check that the parameters are consistent for this sample
        if length(doping_now) ~= length(find(doping_now(1)==doping_now))
            disp(['the doping for sample ' samples{i} ' is not consistent']); 
        end
        if length(thickness_now) ~= length(find(thickness_now(1)==thickness_now))
            disp(['the thickness for sample ' samples{i} ' is not consistent']); 
        end
        [num_meas,columns] = size(raw_now); 
        tau_rev = zeros(num_meas,length(index_FZ)); 
        for k = 1:length(index_FZ)
            SRV_now = SRV_FZ{k};
            for j = 1:num_meas
                [D_now,Dh] = diffusivity(300,'p',doping_now(j),deltan_target);
                W_now = thickness_now(j); 
                try
                    t_index = find(SRV_now(:,1)==raw_now(j,1)); 
                    tau = raw_now(j,2); 
                    tau_surf = (W_now./(2.*SRV_now(t_index,2)))+((1/D_now).*((W_now/pi)^2)); %cm/s
                    tau_rev(j,k) = ((1./tau)-(1./tau_surf))^(-1); 
                catch
                    disp(['There was an error calculating the surface lifetime for ' samples{i} ', time ' num2str(raw_now(j,1)) 's']);
                    tau_rev(j,k) = NaN; 
                end
            end
        end
        FZ_corr{i} = [raw_now(:,1) tau_rev]; 
        if isnan(tau_rev(1,1))==1
            FZ_corr_norm{i} = [raw_now(:,1) tau_rev./tau_rev(1,2)];
        elseif isnan(tau_rev(1,2))==1
            FZ_corr_norm{i} = [raw_now(:,1) tau_rev./tau_rev(1,1)];
        else
            %both were measured?
            FZ_corr_norm{i} = [raw_now(:,1) tau_rev./tau_rev(1,:)];
        end
    end
end

%Now, get the information for the mc-Si as the surface reference. This time
%we won't calculate the SRV but we'll assume similar doping and try
%harmonic sum and ratios. 
index_mc = find(strcmp(samples,'66-2')==1); 
mcSi_harm = cell(size(samples)); 
mcSi_ratio = cell(size(samples)); 
mcSi_harm_norm = cell(size(samples)); 
mcSi_ratio_norm = cell(size(samples)); 
raw_mcSi = lifetime_all{index_mc}; 
for i = 1:length(samples)
    if strcmp(samples{i},'FZ')==0 && strcmp(samples{i},'FZ-new')==0 && ...
            strcmp(samples{i},'FZ-12')==0 && strcmp(samples{i},'66-2')==0 && ...
            strcmp(samples{i},'68-4')==0 && strcmp(samples{i},'60a')==0 && ...
            strcmp(samples{i},'56b')==0 && strcmp(samples{i},'C-1')==0 && ...
            strcmp(samples{i},'H-1')==0 && strcmp(samples{i},'C-2')==0 && ...
            strcmp(samples{i},'H-2')==0
        raw_now = lifetime_all{i};
        [num_meas,columns] = size(raw_now); 
        mcSi_harm_now = zeros(num_meas,1); 
        mcSi_ratio_now = zeros(num_meas,1); 
        for j = 1:num_meas
            t_index = find(raw_mcSi(:,1)==raw_now(j,1)); 
            if isempty(t_index)==1
                disp(['There was no mcSi measurement at time ' num2str(raw_now(j,1)) 's']); 
                mcSi_harm_now(j) = NaN; 
                mcSi_ratio_now(j) = NaN; 
            else
                mcSi_harm_now(j) = 1/((1/raw_now(j,2))-(1/raw_mcSi(t_index,2))); 
                mcSi_ratio_now(j) = raw_now(j,2)/raw_mcSi(t_index,2); 
            end
        end
        mcSi_harm{i} = [raw_now(:,1) mcSi_harm_now]; 
        mcSi_ratio{i} = [raw_now(:,1) mcSi_ratio_now]; 
        mcSi_harm_norm{i} = [raw_now(:,1) mcSi_harm_now./mcSi_harm_now(1)]; 
        mcSi_ratio_norm{i} = [raw_now(:,1) mcSi_ratio_now./mcSi_ratio_now(1)]; 
    end
end

%Now we need to make all of the plots for each plotting group
for i = 1:length(plotting_group)
    group_now = plotting_group{i}; 
    [rows,samp] = size(group_now); 
    lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_FZcorr=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_FZcorr_norm=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_mcSi_harm=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_mcSi_harm_norm=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_mcSi_ratio=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_mcSi_ratio_norm=figure('units','normalized','outerposition',[0 0 1 1]);
    labels_raw = {}; labels_FZcorr = {}; labels_mcSi = {}; 
    hFZ = []; hFZnorm = [];
    cm = colormap(hsv(samp));
    for j = 1:samp
        index = find(strcmp(group_now{1,j},samples)==1);
        raw_now = lifetime_all{index}; 
        norm_now = norm_lifetime_all{index};
        FZ_corr_now = FZ_corr{index}; 
        FZ_corr_norm_now = FZ_corr_norm{index}; 
        mcSi_harm_now = mcSi_harm{index}; 
        mcSi_harm_norm_now = mcSi_harm_norm{index}; 
        mcSi_ratio_now = mcSi_ratio{index}; 
        mcSi_ratio_norm_now = mcSi_ratio_norm{index}; 
        figure(lifetime_raw); 
        plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
        hold all; 
        figure(lifetime_norm); 
        plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
        hold all;  
        labels_raw{end+1} = group_now{2,j};  
        if strcmp(group_now{1,j},'FZ')==0 && strcmp(group_now{1,j},'FZ-new')==0 && ...
            strcmp(group_now{1,j},'68-4')==0 && strcmp(group_now{1,j},'60a')==0 && ...
            strcmp(group_now{1,j},'56b')==0
            for k = 1:length(index_FZ)
                figure(lifetime_FZcorr); 
                if k == 1
                    hFZ(end+1)=plot(FZ_corr_now(:,1),FZ_corr_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:)); 
                    hold on; 
                    figure(lifetime_FZcorr_norm); 
                    hFZnorm(end+1)=plot(FZ_corr_norm_now(:,1),FZ_corr_norm_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:));
                    hold on;  
                else
                    plot(FZ_corr_now(:,1),FZ_corr_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:)); 
                    hold on; 
                    figure(lifetime_FZcorr_norm); 
                    plot(FZ_corr_norm_now(:,1),FZ_corr_norm_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:));
                    hold on;  
                end
            end
            labels_FZcorr{end+1} = group_now{2,j}; 
        end
        if strcmp(group_now{1,j},'FZ')==0 && strcmp(group_now{1,j},'FZ-new')==0 && ...
            strcmp(group_now{1,j},'FZ-12')==0 && strcmp(group_now{1,j},'66-2')==0 && ...
            strcmp(group_now{1,j},'68-4')==0 && strcmp(group_now{1,j},'60a')==0 && ...
            strcmp(group_now{1,j},'56b')==0 && strcmp(group_now{1,j},'C-1')==0 && ...
            strcmp(group_now{1,j},'C-2')==0 && strcmp(group_now{1,j},'H-1')==0 && ...
            strcmp(group_now{1,j},'H-2')==0
            figure(lifetime_mcSi_harm); 
            plot(mcSi_harm_now(:,1),mcSi_harm_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
            hold all; 
            figure(lifetime_mcSi_harm_norm); 
            plot(mcSi_harm_norm_now(:,1),mcSi_harm_norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
            hold all;
            figure(lifetime_mcSi_ratio); 
            plot(mcSi_ratio_now(:,1),mcSi_ratio_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
            hold all; 
            figure(lifetime_mcSi_ratio_norm); 
            plot(mcSi_ratio_norm_now(:,1),mcSi_ratio_norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
            hold all;
            labels_mcSi{end+1} = group_now{2,j}; 
        end
    end
    figure(lifetime_raw); 
    xlabel('time [s]','FontSize',25); 
    ylabel('lifetime [s]','FontSize',25); 
    legend(labels_raw); 
    title(plotting_names{i},'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_raw,[savedirname '\' plotting_names{i} '_raw' savename]);
    print(lifetime_raw,'-dpng','-r0',[savedirname '\' plotting_names{i} '_raw' savename '.png']);
    
    figure(lifetime_norm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('normalized lifetime [-]','FontSize',25); 
    axis([0 max_time 0 2]);
    legend(labels_raw); 
    title(plotting_names{i},'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_norm,[savedirname '\' plotting_names{i} '_norm' savename]);
    print(lifetime_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_norm' savename '.png']);
    
    figure(lifetime_FZcorr); 
    xlabel('time [s]','FontSize',25); 
    ylabel('lifetime [s]','FontSize',25); 
    legend(hFZ,labels_FZcorr); 
    title([plotting_names{i} ' corrected by FZ SRV'],'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_FZcorr,[savedirname '\' plotting_names{i} '_FZcorr' savename]);
    print(lifetime_FZcorr,'-dpng','-r0',[savedirname '\' plotting_names{i} '_FZcorr' savename '.png']);
    
    figure(lifetime_FZcorr_norm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('normalized lifetime [-]','FontSize',25); 
    axis([0 max_time 0 2]);
    legend(hFZnorm,labels_FZcorr); 
    title([plotting_names{i} ' corrected by FZ SRV'],'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_FZcorr_norm,[savedirname '\' plotting_names{i} '_FZcorr_norm' savename]);
    print(lifetime_FZcorr_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_FZcorr_norm' savename '.png']);
    
    figure(lifetime_mcSi_harm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('lifetime [s]','FontSize',25); 
    legend(labels_mcSi); 
    title([plotting_names{i} ' corrected by mcSi control harmonic sum'],'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_mcSi_harm,[savedirname '\' plotting_names{i} '_mcSi_harm' savename]);
    print(lifetime_mcSi_harm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_harm' savename '.png']);
    
    figure(lifetime_mcSi_harm_norm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('normalized lifetime [-]','FontSize',25); 
    axis([0 max_time 0 2]);
    legend(labels_mcSi); 
    title([plotting_names{i} ' corrected by mcSi control harmonic sum'],'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_mcSi_harm_norm,[savedirname '\' plotting_names{i} '_mcSi_harm_norm' savename]);
    print(lifetime_mcSi_harm_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_harm_norm' savename '.png']);
    
    figure(lifetime_mcSi_ratio); 
    xlabel('time [s]','FontSize',25); 
    ylabel('lifetime [s]','FontSize',25); 
    legend(labels_mcSi); 
    title([plotting_names{i} ' corrected by mcSi control ratio'],'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_mcSi_ratio,[savedirname '\' plotting_names{i} '_mcSi_ratio' savename]);
    print(lifetime_mcSi_ratio,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_ratio' savename '.png']);
    
    figure(lifetime_mcSi_ratio_norm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('normalized lifetime [-]','FontSize',25); 
    axis([0 max_time 0 2]);
    legend(labels_mcSi); 
    title([plotting_names{i} ' corrected by mcSi control ratio'],'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_mcSi_ratio_norm,[savedirname '\' plotting_names{i} '_mcSi_ratio_norm' savename]);
    print(lifetime_mcSi_ratio_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_ratio_norm' savename '.png']);
end

%Let's save all of the processed data as a .mat file for easy use later
save([savedirname '\' bora '_all_data' savename '.mat'],'lifetime_all',...
    'norm_lifetime_all','FZ_corr','FZ_corr_norm','mcSi_harm',...
    'mcSi_harm_norm','mcSi_ratio','mcSi_ratio_norm','doping_all',...
    'thickness_all','samples','SRV_FZ','surface_control');

