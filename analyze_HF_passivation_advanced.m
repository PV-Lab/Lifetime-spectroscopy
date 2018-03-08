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
bora = 'set-a'; %'set-b' or 'set-a' or 'compE' or 'compare' if you want to compare sets a and b directly
%Most recent directory that we want to analyze now. 
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 9 2017';
%where we want to save any new, non-sample-specific data
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\PVSC abstract\lifetime to plot'; 
%Spreadsheet specification for the actual measurements
spreadsheet = 'new'; %old (before TS) or new (after TS)
%target injection level for the measurements, used to make degradation
%curves
deltan_target = 6e14; %cm^-3
%Error in lifetime measurement (approximate based on 66-2, valid for mc-si
%samples only)
lifetime_error = 0.11;
type = 'p'; 

%-----------------------------
%which sample set, b or a? 
[samples,dirnames,labels,savename,surface_control,...
    plotting_group,plotting_names,meas_details,max_time,...
    latesta,latestb,lifetime_analysis]=HFbora(bora);

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
                deltan = datanow(:,1); tau = datanow(:,2); 
                cutoff = 5e15; 
                [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
                %We might always want to remove some low injection data
%                 disp('Select the region for cutting off the LOW injection data');
%                 [cutoff,nothing]=ginput(1);
                cutoff = 1e14; 
                [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
                curves(count)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2,'color',cm(k,:)); 
                hold all; 
                label{count} = labels{k};
                xlim([5e13 1e17])
                count = count+1; 
                lifetime_store{k,i} = [deltan_rev,tau_rev];
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
% 9) Nt*, calculated from the harmonic difference between each subsequent
% measurement and the initial measurement. (Ntstar)

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
Ntstar = cell(length(samples),1); 
for i = 1:length(samples)
    meas_thissample = meas(i,:);
    lifetime_store = [];
    doping_store = [];
    thickness_store = []; 
    Ntstar_thissample = []; 
    Ntstar_error = []; 
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
            %calculate and store the harmonic difference = Ntstar
            Ntstar_thissample(j,1) = (1/lifetime_store(j,1))-(1/lifetime_store(1,1)); %this should be in units of inverse seconds
            %Calculate the error associated with Ntstar
            dNtdtdeg = 1/(lifetime_store(j,1)^2); 
            dNtdtinit = 1/(lifetime_store(1,1)^2);
            dtdeg = lifetime_error*lifetime_store(j,1); 
            dtinit = lifetime_error*lifetime_store(j,1); 
            Ntstar_error(j,1) = sqrt(((dNtdtdeg*dtdeg)^2)+((dNtdtinit*dtinit)^2)); 
        end
    end
    nan_indices = find(isnan(meas_thissample)==1); 
    meas_thissample(nan_indices) = []; 
    lifetime_all{i} = [meas_thissample' lifetime_store]; 
    norm_lifetime_all{i} = [meas_thissample' lifetime_store./lifetime_store(1)];
    thickness_all{i} = thickness_store; 
    doping_all{i} = doping_store; 
    Ntstar{i} = [meas_thissample' Ntstar_thissample Ntstar_error]; 
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
            strcmp(samples{i},'56b')==0 && strcmp(samples{i},'FZ-new2')==0
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
            strcmp(samples{i},'H-2')==0 && strcmp(samples{i},'FZ-new2')==0
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
    Ntstar_time=figure('units','normalized','outerposition',[0 0 1 1]);
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
        Ntstar_now = Ntstar{index}; 
        figure(Ntstar_time); 
        plot(Ntstar_now(:,1),Ntstar_now(:,2),'o','LineWidth',3,'MarkerSize',10);
        hold all; 
        figure(lifetime_raw); 
        plot(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
        hold all; 
        figure(lifetime_norm); 
        plot(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
        hold all;  
        labels_raw{end+1} = group_now{2,j};  
        if strcmp(group_now{1,j},'FZ')==0 && strcmp(group_now{1,j},'FZ-new')==0 && ...
            strcmp(group_now{1,j},'68-4')==0 && strcmp(group_now{1,j},'60a')==0 && ...
            strcmp(group_now{1,j},'56b')==0 && strcmp(group_now{1,j},'FZ-new2')==0
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
            strcmp(group_now{1,j},'H-2')==0 && strcmp(group_now{1,j},'FZ-new2')==0
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
    
    figure(Ntstar_time); 
    xlabel('time [s]','FontSize',25); 
    ylabel('N_t^* [s^-^1]','FontSize',25); 
    legend(labels_raw); 
    title(plotting_names{i},'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(Ntstar_time,[savedirname '\' plotting_names{i} '_Ntstar' savename]);
    print(Ntstar_time,'-dpng','-r0',[savedirname '\' plotting_names{i} '_Ntstar' savename '.png']);
    
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
    'thickness_all','samples','SRV_FZ','surface_control','Ntstar');

%% Make the degradation curves comparing data already processed - compare sets a and b
close all; clc; 

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
        %load the right data
        if strcmp(group_now{2,j},'a')==1
            load(latesta); 
        elseif strcmp(group_now{2,j},'b')==1
            load(latestb); 
        else
            disp('error loading data'); 
        end
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
            strcmp(group_now{1,j},'56b')==0 && strcmp(group_now{1,j},'FZ-new2')==0
            [num_meas,index_FZ] = size(FZ_corr_now); 
            index_FZ = index_FZ-1; %the first column is the time
            for k = 1:index_FZ
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
            strcmp(group_now{1,j},'H-2')==0 && strcmp(group_now{1,j},'FZ-new2')==0
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
    xlim([0 max_time]);
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
    xlim([0 max_time]);
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
    xlim([0 max_time]);
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
    xlim([0 max_time]);
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
%% Lifetime spectroscopy analysis
%Currently the inputs for this are only set for set-a
close all; clc; 
fit_tries = 1e5;
%For each sample, we will have a set of measurements - we want to fit each
%one 
[meas,samples] = xlsread(meas_details,'measurements');
samples(1,:) = []; 
[times,filenames] = xlsread(meas_details,'filenames'); 
%Make these the same size
filenames = filenames(2:end,2); 

%Figure out the number of samples we are going to analyze this round 
[rows,num_samples] = size(lifetime_analysis); 

%Initialize the structures
lifetime_all = cell(num_samples,2); 
thickness_all = cell(num_samples,1); 
doping_all = cell(num_samples,1); 

%loop over only the samples we want for lifetime analysis
for i = 1:length(samples)
    %figure out which measurements were taking for this sample
%     index = find(strcmp(samples,lifetime_analysis{1,i})==1); 
    index = i; 
    meas_thissample = meas(index,:); 
    lifetime_store = {};
    doping_store = [];
    thickness_store = [];
    %We will do the same procedure for each measurement at each degradation
    %time
    for j = 1:length(meas_thissample)
        if isnan(meas_thissample(j))==0
            %now create the proper filename
            findex = find(meas_thissample(j)==times);  
            filename = [filenames{findex} '\' samples{index} '\Raw_data.mat'];
            load(filename);
            filename = [filenames{findex} '\' samples{index} '\meas_info.mat'];
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
            lifetime_store{j,1} = [deltan,tau]; 
        end
    end
    thickness_all{i} = thickness_store; 
    doping_all{i} = doping_store; 
    nan_indices = find(isnan(meas_thissample)==1); 
    meas_thissample(nan_indices) = []; 
    lifetime_all{i,1} = meas_thissample';
    lifetime_all{i,2} = lifetime_store;
end
            
%Now we need to make the corrections - first, using FZ as the surface
%reference.
index_FZ = [];
for i = 1:length(surface_control)
    index_FZ(end+1) = find(strcmp(samples,surface_control{i})==1); 
end
%We need to get the surface information from FZ before we loop
%This needs to be done in an injection dependent way
SRV_FZ = cell(length(index_FZ),2); 
SRV_label = cell(size(index_FZ)); 
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
    raw_FZ = lifetime_all{index_FZ(i),2};
    meas = lifetime_all{index_FZ(i),1}; 
    SRV_thisFZ = cell(length(meas),1); 
    for j = 1:length(meas)
        raw_FZ_data = raw_FZ{j};
        [injection,col] = size(raw_FZ_data); 
        SRV = zeros(injection,1); 
        for k = 1:injection
            [De,Dh] = diffusivity(300,'p',doping_FZ(j),raw_FZ_data(k,1)); 
            D_FZ  = De; %cm2/s
            tau_intr = Richter(300,raw_FZ_data(k,1),doping_FZ(j),'p');
            tau_surf = ((1./raw_FZ_data(k,2))-(1./tau_intr)).^(-1);
            SRV(k) = thickness_FZ(j)./((tau_surf-((1/D_FZ)*((thickness_FZ(j)/pi)^2))).*2);
        end
        %get rid of any negative or nan SRV's
        index_neg = find(SRV<0); 
        deltan = raw_FZ_data(:,1); 
        SRV(index_neg) = []; deltan(index_neg) = []; 
        index_nan = find(isnan(SRV)==1); 
        SRV(index_nan) = []; deltan(index_nan) = []; 
        SRV_thisFZ{j,1} = [deltan SRV]; 
    end
    SRV_FZ{i,1} = meas;
    SRV_FZ{i,2} = SRV_thisFZ; 
    SRV_label{i} = samples{index_FZ(i)}; 
end
    
%Now we have the lifetime for every sample, and the FZ lifetime for every
%time. We need to bring them together to get the injection dependent SRH
%lifetime, which we can then fit, for each sample of interest at each time.
for i = 1:num_samples
    if strcmp(lifetime_analysis{1,i},'FZ')==0 && strcmp(lifetime_analysis{1,i},'FZ-new')==0 && ...
            strcmp(lifetime_analysis{1,i},'68-4')==0 && strcmp(lifetime_analysis{1,i},'60a')==0 && ...
            strcmp(lifetime_analysis{1,i},'56b')==0 && strcmp(lifetime_analysis{1,i},'FZ-new2')==0
        %Then we explicitly find the lifetime using the FZ wafer as a
        %reference
        %....
    elseif strcmp(lifetime_analysis{1,i},'68-4')==1 || ...
            strcmp(lifetime_analysis{1,i},'60a')==1 || ...
            strcmp(lifetime_analysis{1,i},'56b')==1
        %Then we will go ahead and use the harmonic sum 
        %....
    end
end

    
    
    
    
    
%% 
    
    lifetimefig=figure; 
    %Load the initial lifetime, which should be the first entry in the
    %dirnames list
    try 
        initial_data=load([dirnames{1} '\' lifetime_analysis{1,i} '\Raw_data.mat']); 
        load([dirnames{1} '\' lifetime_analysis{1,i} '\meas_info.mat']); 
        initial_info = info; 
    catch
        %If that's not successful, just tell me that the sample doesn't
        %exist for that dirname. There will be a script error. 
        warning(['Error accessing data for directory ' num2str(1) ', sample ' lifetime_analysis{1,i}]);
    end
    if length(initial_data.dataSave)>1
        if length(initial_data.dataSave)>2
            t = 3; 
        else
            t = 2; 
        end
    else
        t = 1; 
    end
    data = initial_data.dataSave{t}; 
    deltan_initial=data(:,1);tau_initial=data(:,2);
    figure(lifetimefig); loglog(deltan_initial,tau_initial); 
    %Now we have: initial_data.fileListShort, initial_data.dataSave
    %initial_info.thickness,resistivity,doping,temperature
    %Now we need to loop over each measured lifetime and calculate the SRH
    %lifetime by the harmonic difference
    for j = 2:length(dirnames)
        %Load the data for this sample
        try 
            load([dirnames{j} '\' lifetime_analysis{1,i} '\Raw_data.mat']); 
            load([dirnames{j} '\' lifetime_analysis{1,i} '\meas_info.mat']); 
            flag = 1; 
        catch
            %If that's not successful, just tell me that the sample doesn't
            %exist for that dirname. There will be a script error. 
            warning(['Error accessing data for directory ' num2str(j) ', sample ' lifetime_analysis{1,i}]);
            flag = 0; 
        end
        %If we've successfully loaded the data, we need to do some
        %operations with it
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
            for k  = t %just the second measurement in each set
                data = dataSave{k}; 
                deltan = data(:,1); tau = data(:,2); 
                figure(lifetimefig); hold all; loglog(deltan,tau); 
                %We need to interpolate this data set at the injection
                %levels for our other sample
                try
                    [tau] = interp1(deltan,tau,deltan_initial); 
                catch
                    [deltan,tau] = remove_duplicates(deltan,tau);
                    try
                        [tau] = interp1(deltan,tau,deltan_initial); 
                    catch
                        [deltan,tau] = remove_duplicates(deltan,tau);
                        [tau] = interp1(deltan,tau,deltan_initial); 
                    end
                end
                %If we've exited successfully deltan should now be
                %deltan_initial
                deltan = deltan_initial; 
                %Ask where to crop the measured lifetime
                figure;
                loglog(deltan,tau,'-o');
%                 disp('Select the region for cutting off the HIGH injection data');
%                 [cutoff,nothing]=ginput(1);
                cutoff = 5e15; 
                [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
                [deltan_initial_rev,tau_initial_rev] = remove_highinj(deltan_initial,tau_initial,cutoff); 
                %We might always want to remove some low injection data
%                 disp('Select the region for cutting off the LOW injection data');
%                 [cutoff,nothing]=ginput(1);
                cutoff = 1e14; 
                [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
                [deltan_initial_rev,tau_initial_rev] = remove_lowinj(deltan_initial_rev,tau_initial_rev,cutoff); 
                hold all;
                loglog(deltan_rev,tau_rev,'+');
                legend('Before cutoff','After cutoff'); 
                %Now we calculate the harmonic difference, which represents
                %the SRH lifetime
                tauSRH = ((1./tau_rev)-(1./tau_initial_rev)).^(-1); %units of seconds
                %Linearize the lifetime
                %Get sample parameters at specified temperature
                [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(info(t).temperature+273.15,info(t).doping,type); 
                %Normalized carrier density
                if type == 'p'
                    X = (n0+deltan_rev)./(p0+deltan_rev);
                elseif type == 'n'
                    X = (p0+deltan_rev)./(n0+deltan_rev);
                end
                %Send the lifetime for fitting
                indices = find(tauSRH<0 | isnan(tauSRH)==1); 
                tauSRH(indices) = [];
                X(indices) = [];
                xlswrite([savedirname '\Linearized_data_' lifetime_analysis{1,i} '.xlsx'],[X,tauSRH],['Sheet' num2str(j)]); 
                [one_defect{j,i},MSE_one{j,i},two_defects{j,i},MSE_two{j,i},three_defects{j,i},MSE_three{j,i},all_parameters_store{j,i},all_MSE_store{j,i}] = fit_murphy_master(X,tauSRH.*1e6,info(t).temperature,savedirname,fit_tries);
                to_write = zeros(6,3); 
                to_write(1:2,1) = one_defect{j,i}';
                twodef = two_defects{j,i}; 
                to_write(1:2,2) = twodef(1,:)';
                to_write(3:4,2) = twodef(2,:)';
                threedef = three_defects{j,i}; 
                to_write(1:2,3) = threedef(1,:)';
                to_write(3:4,3) =threedef(2,:)';
                to_write(5:6,3) = threedef(3,:)';
                xlswrite([savedirname '\Linearized_data_' lifetime_analysis{1,i} '.xlsx'],to_write,['Sheet' num2str(j)],'C1:E6'); 
                X_store{j,i} = X; 
                tau_SRH_store{j,i} = tauSRH; 
            end
        end
    end
    figure(lifetimefig); 
    xlabel('excess carrier density [cm^-^3]'); 
    ylabel('lifetime [s]'); 
    legend(labels'); 
    title([lifetime_analysis{1,i} ': ' lifetime_analysis{2,i}]); 
    hgsave(lifetimefig,[savedirname '\' lifetime_analysis{1,i} '_meas_lifetimes']);
    print(lifetimefig,'-dpng','-r0',[savedirname '\' lifetime_analysis{1,i} '_meas_lifetimes.png']); 
end
best_fits = struct('one_defect',one_defect,'MSE_one',MSE_one,...
    'two_defects',two_defects,'MSE_two',MSE_two,'three_defects',...
    three_defects,'MSE_three',MSE_three,'all_fits',all_parameters_store,...
    'all_MSE',all_MSE_store);
save([savedirname '\best_fits.mat'],'best_fits','lifetime_analysis');
save([savedirname '\original_linearized.mat'],'X_store','tau_SRH_store','lifetime_analysis');

%Outside of this, need to modify the fit parameters in Excel, match
%defects. Then need to bring these into MATLAB and calculate the E-k
%values, extract midgap value, and then extract for plotting in Origin. 

%% Making E-k curves - from Excel
close all; clc; 

for i = 1:length(lifetime_analysis)
    %Load the first measurement just so we can get the relevant parameters
    load([dirnames{1} '\' lifetime_analysis{1,i} '\meas_info.mat']); 
    %Get the data from Excel
    filename = [savedirname '\IDLS two and three curve fitting_' lifetime_analysis{1,i} '.xlsm'];
    fits = xlsread(filename,'Summary','B3:E20'); %Note - they should be sorted in increasing time
    [m,n] = size(fits); 
    %Make some figures for each defect
    def1=figure; def2=figure;
    %for each time for this sample 
    for j = 1:m
        %Get the fit at this time
        defect_1 = fits(j,1:2); 
        defect_2 = fits(j,3:4); 
        %Make the E-k curves based on these new fits
        [Et{j,1},k{j,1},alphanN{j,1}]=generate_Ek(defect_1,info(1).temperature+273.15,info(1).doping,type);
        [Et{j,2},k{j,2},alphanN{j,2}]=generate_Ek(defect_2,info(1).temperature+273.15,info(1).doping,type);
        %Plot the E-k curves for both defects
        figure(def1);
        plot(Et{j,1},k{j,1},'-','LineWidth',2); 
        hold all;
        figure(def2); 
        plot(Et{j,2},k{j,2},'-','LineWidth',2); 
        hold all;
        %Get the value at the intrinsic point 
        index = find(abs(Et{j,1}-0)==min(abs(Et{j,1}-0))); 
        k_intrinsic_store(j,1) = k{j,1}(index);
        alphanN_intrinsic_store(j,1) = alphanN{j,1}(index); 
        index = find(abs(Et{j,2}-0)==min(abs(Et{j,2}-0))); 
        k_intrinsic_store(j,2) = k{j,2}(index);
        alphanN_intrinsic_store(j,2) = alphanN{j,2}(index); 
    end
    %format the figurse and save
    figure(def1); 
    xlabel('E_t-E_i [eV]','FontSize',20); 
    ylabel('k [-]','FontSize',20);
    legend(labels'); 
    title([lifetime_analysis{1,i} ': Defect 1'],'FontSize',30); 
    hgsave(def1,[savedirname '\' lifetime_analysis{1,i} '_defect1_E-k']);
    print(def1,'-dpng','-r0',[savedirname '\' lifetime_analysis{1,i} '_defect1_E-k.png']); 
    figure(def2); 
    xlabel('E_t-E_i [eV]','FontSize',20); 
    ylabel('k [-]','FontSize',20);
    legend(labels(2:end)'); 
    title([lifetime_analysis{1,i} ': Defect 2'],'FontSize',30); 
    hgsave(def1,[savedirname '\' lifetime_analysis{1,i} '_defect2_E-k']);
    print(def1,'-dpng','-r0',[savedirname '\' lifetime_analysis{1,i} '_defect2_E-k.png']); 
    %Store the data
    all_kintr{i} = k_intrinsic_store; 
    all_alphanNintr{i} = alphanN_intrinsic_store; 
    Et_store{i} = Et; 
    k_store{i} = k; 
    alphanN_store{i} = alphanN; 
end

%% Making E-k curves - directly from .mat
close all; clc; 
load([savedirname '\best_fits.mat']); 
% cutoff_X = 0.025;
[rows,num_samples] = size(lifetime_analysis); 
%loop over only the samples we want for lifetime analysis
for i = 1:num_samples
    %Load the first measurement just so we can get the relevant parameters
    load([dirnames{1} '\' lifetime_analysis{1,i} '\meas_info.mat']); 
%     %Get the data from Excel
%     fits = best_fits; 
%     [m,n] = size(fits); 
    %Make some figures for each defect
    def1=figure; def2=figure;
    %for each time for this sample 
    for j = 2:length(dirnames)
%         index = find(abs(X{j}-cutoff_X)==min(abs(X{j}-cutoff_X))); 
%         target_tau = tau_SRH_store{j}(index); 
        fits = best_fits(j).two_defects;
        %Get the fit at this time
        if fits(1,1)<fits(2,1)
            defect_1 = fits(2,:);
            defect_2 = fits(1,:);
        else
            defect_1 = fits(1,:);
            defect_2 = fits(2,:); 
        end
        %Make the E-k curves based on these new fits
        [Et{j,1},k{j,1},alphanN{j,1}]=generate_Ek(defect_1,info(1).temperature+273.15,info(1).doping,type);
        [Et{j,2},k{j,2},alphanN{j,2}]=generate_Ek(defect_2,info(1).temperature+273.15,info(1).doping,type);
        %Plot the E-k curves for both defects
        figure(def1);
        plot(Et{j,1},k{j,1},'-','LineWidth',2); 
        hold all;
        figure(def2); 
        plot(Et{j,2},k{j,2},'-','LineWidth',2); 
        hold all;
        %Get the value at the intrinsic point 
        index = find(abs(Et{j,1}-0)==min(abs(Et{j,1}-0))); 
        k_intrinsic_store(j,1) = k{j,1}(index);
        alphanN_intrinsic_store(j,1) = alphanN{j,1}(index); 
        index = find(abs(Et{j,2}-0)==min(abs(Et{j,2}-0))); 
        k_intrinsic_store(j,2) = k{j,2}(index);
        alphanN_intrinsic_store(j,2) = alphanN{j,2}(index); 
    end
    %format the figurse and save
    figure(def1); 
    xlabel('E_t-E_i [eV]','FontSize',20); 
    ylabel('k [-]','FontSize',20);
    legend(labels(2:end)'); 
    title([lifetime_analysis{1,i} ': Defect 1'],'FontSize',30); 
    hgsave(def1,[savedirname '\' lifetime_analysis{1,i} '_defect1_E-k']);
    print(def1,'-dpng','-r0',[savedirname '\' lifetime_analysis{1,i} '_defect1_E-k.png']); 
    figure(def2); 
    xlabel('E_t-E_i [eV]','FontSize',20); 
    ylabel('k [-]','FontSize',20);
    legend(labels'); 
    title([lifetime_analysis{1,i} ': Defect 2'],'FontSize',30); 
    hgsave(def1,[savedirname '\' lifetime_analysis{1,i} '_defect2_E-k']);
    print(def1,'-dpng','-r0',[savedirname '\' lifetime_analysis{1,i} '_defect2_E-k.png']); 
    %Store the data
    all_kintr{i} = k_intrinsic_store; 
    all_alphanNintr{i} = alphanN_intrinsic_store; 
    Et_store{i} = Et; 
    k_store{i} = k; 
    alphanN_store{i} = alphanN; 
end
    