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
bora = 'follow-up'; %'set-b' or 'set-a' or 'compE' or 'compare' if you want to compare sets a and b directly
%Most recent directory that we want to analyze now. 
% dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set b\113020s';
% dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set a\3761830s';
% dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\compE\1518840s';

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\Follow-up experiment\lifetime data\March 2018 as received UKN\degradation\May 10 2018\457000s';

%where we want to save any new, non-sample-specific data
% savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set b\113020s\lifetime spectroscopy'; 
% savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set a\3761830s\lifetime spectroscopy';
% savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\compE\1518840s\lifetime spectroscopy';
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\Follow-up experiment\lifetime data\March 2018 as received UKN\degradation\summary\464974s\1-64'; 
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
                    t = 1; %1-64
%                     t = 3; %1-1 high
                else
                    t = 2; %high measurement
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
                if length(dataSave)>2
                    t = 1; %1-64
%                     t = 3; %1-1 high
                else
                    t = 2; %high measurement
                end
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
    if strcmp(samples{i},'59b')==1
        meas_thissample=meas_thissample+((6*60*60)+8*60);
    end
    nan_indices = find(isnan(meas_thissample)==1); 
    meas_thissample(nan_indices) = []; 
    lifetime_all{i} = [meas_thissample' lifetime_store]; 
    norm_lifetime_all{i} = [meas_thissample' lifetime_store./lifetime_store(1)];
    thickness_all{i} = thickness_store; 
    doping_all{i} = doping_store; 
    Ntstar{i} = [meas_thissample' Ntstar_thissample Ntstar_error]; 
end

% %Now we need to make the corrections - first, using FZ as the surface
% %reference.
% index_FZ = [];
% for i = 1:length(surface_control)
%     index_FZ(end+1) = find(strcmp(samples,surface_control{i})==1); 
% end
% %We need to get the surface information from FZ before we loop
% SRV_FZ = cell(size(index_FZ)); 
% SRV_label = cell(size(index_FZ)); 
% figure('units','normalized','outerposition',[0 0 1 1]);; 
% for i = 1:length(index_FZ)
%     doping_FZ = doping_all{index_FZ(i)}; %cm-3
%     thickness_FZ = thickness_all{index_FZ(i)}; %cm
%     %Check that the parameters are consistent for this sample
%     if length(doping_FZ) ~= length(find(doping_FZ(1)==doping_FZ))
%         disp(['the doping for sample ' samples{index_FZ(i)} ' is not consistent']); 
%     end
%     if length(thickness_FZ) ~= length(find(thickness_FZ(1)==thickness_FZ))
%         disp(['the thickness for sample ' samples{index_FZ(i)} ' is not consistent']); 
%     end
%     raw_FZ = lifetime_all{index_FZ(i)}; 
%     [meas,nothing]=size(raw_FZ); 
%     SRV = zeros(meas,1); 
%     for j = 1:length(SRV)
%         [De,Dh] = diffusivity(300,'p',doping_FZ(j),deltan_target); 
%         D_FZ  = De; %cm2/s
%         tau_intr = Richter(300,deltan_target,doping_FZ(j),'p');
%         tau_surf = ((1./raw_FZ(j,2))-(1./tau_intr)).^(-1);
%         SRV(j) = thickness_FZ(j)./((tau_surf-((1/D_FZ)*((thickness_FZ(j)/pi)^2))).*2);
%     end
%     SRV_FZ{i} = [raw_FZ(:,1) SRV]; 
%     if raw_FZ(1,1) == 0
%         raw_FZ(1,1) = 1; 
%     end
%     semilogx(raw_FZ(:,1),SRV,'LineWidth',3); 
%     hold all; 
%     SRV_label{i} = samples{index_FZ(i)}; 
% end
% axis([1 max_time 0 10]);
% xlabel('time [s]','FontSize',25); 
% ylabel('SRV [cm/s]','FontSize',25); 
% legend(SRV_label); 
% set(0,'defaultAxesFontSize', 20)
% hgsave(gcf,[savedirname '\SRV' savename]);
% print(gcf,'-dpng','-r0',[savedirname '\SRV' savename '.png']);
% 
% FZ_corr = cell(size(samples)); 
% FZ_corr_norm = cell(size(samples)); 
% count = 1; 
% for i = 1:length(samples)
%     if strcmp(samples{i},'FZ')==0 && strcmp(samples{i},'FZ-new')==0 && ...
%             strcmp(samples{i},'68-4')==0 && strcmp(samples{i},'60a')==0 && ...
%             strcmp(samples{i},'56b')==0 && strcmp(samples{i},'FZ-new2')==0
%         raw_now = lifetime_all{i};
%         doping_now = doping_all{i}; 
%         thickness_now = thickness_all{i}; 
%         %Check that the parameters are consistent for this sample
%         if length(doping_now) ~= length(find(doping_now(1)==doping_now))
%             disp(['the doping for sample ' samples{i} ' is not consistent']); 
%         end
%         if length(thickness_now) ~= length(find(thickness_now(1)==thickness_now))
%             disp(['the thickness for sample ' samples{i} ' is not consistent']); 
%         end
%         [num_meas,columns] = size(raw_now); 
%         tau_rev = zeros(num_meas,length(index_FZ)); 
%         for k = 1:length(index_FZ)
%             SRV_now = SRV_FZ{k};
%             for j = 1:num_meas
%                 [D_now,Dh] = diffusivity(300,'p',doping_now(j),deltan_target);
%                 W_now = thickness_now(j); 
%                 try
%                     t_index = find(SRV_now(:,1)==raw_now(j,1)); 
%                     tau = raw_now(j,2); 
%                     tau_surf = (W_now./(2.*SRV_now(t_index,2)))+((1/D_now).*((W_now/pi)^2)); %cm/s
%                     tau_rev(j,k) = ((1./tau)-(1./tau_surf))^(-1); 
%                 catch
%                     disp(['There was an error calculating the surface lifetime for ' samples{i} ', time ' num2str(raw_now(j,1)) 's']);
%                     tau_rev(j,k) = NaN; 
%                 end
%             end
%         end
%         FZ_corr{i} = [raw_now(:,1) tau_rev]; 
%         if isnan(tau_rev(1,1))==1
%             FZ_corr_norm{i} = [raw_now(:,1) tau_rev./tau_rev(1,2)];
%         elseif isnan(tau_rev(1,2))==1
%             FZ_corr_norm{i} = [raw_now(:,1) tau_rev./tau_rev(1,1)];
%         else
%             %both were measured?
%             FZ_corr_norm{i} = [raw_now(:,1) tau_rev./tau_rev(1,:)];
%         end
%     end
% end
% 
% %Now, get the information for the mc-Si as the surface reference. This time
% %we won't calculate the SRV but we'll assume similar doping and try
% %harmonic sum and ratios. 
% index_mc = find(strcmp(samples,'66-2')==1); 
% mcSi_harm = cell(size(samples)); 
% mcSi_ratio = cell(size(samples)); 
% mcSi_harm_norm = cell(size(samples)); 
% mcSi_ratio_norm = cell(size(samples)); 
% raw_mcSi = lifetime_all{index_mc}; 
% for i = 1:length(samples)
%     if strcmp(samples{i},'FZ')==0 && strcmp(samples{i},'FZ-new')==0 && ...
%             strcmp(samples{i},'FZ-12')==0 && strcmp(samples{i},'66-2')==0 && ...
%             strcmp(samples{i},'68-4')==0 && strcmp(samples{i},'60a')==0 && ...
%             strcmp(samples{i},'56b')==0 && strcmp(samples{i},'C-1')==0 && ...
%             strcmp(samples{i},'H-1')==0 && strcmp(samples{i},'C-2')==0 && ...
%             strcmp(samples{i},'H-2')==0 && strcmp(samples{i},'FZ-new2')==0
%         raw_now = lifetime_all{i};
%         [num_meas,columns] = size(raw_now); 
%         mcSi_harm_now = zeros(num_meas,1); 
%         mcSi_ratio_now = zeros(num_meas,1); 
%         for j = 1:num_meas
%             t_index = find(raw_mcSi(:,1)==raw_now(j,1)); 
%             if isempty(t_index)==1
%                 disp(['There was no mcSi measurement at time ' num2str(raw_now(j,1)) 's']); 
%                 mcSi_harm_now(j) = NaN; 
%                 mcSi_ratio_now(j) = NaN; 
%             else
%                 mcSi_harm_now(j) = 1/((1/raw_now(j,2))-(1/raw_mcSi(t_index,2))); 
%                 mcSi_ratio_now(j) = raw_now(j,2)/raw_mcSi(t_index,2); 
%             end
%         end
%         mcSi_harm{i} = [raw_now(:,1) mcSi_harm_now]; 
%         mcSi_ratio{i} = [raw_now(:,1) mcSi_ratio_now]; 
%         mcSi_harm_norm{i} = [raw_now(:,1) mcSi_harm_now./mcSi_harm_now(1)]; 
%         mcSi_ratio_norm{i} = [raw_now(:,1) mcSi_ratio_now./mcSi_ratio_now(1)]; 
%     end
% end

%Now we need to make all of the plots for each plotting group
for i = 1:length(plotting_group)
    group_now = plotting_group{i}; 
    [rows,samp] = size(group_now); 
    lifetime_raw=figure('units','normalized','outerposition',[0 0 1 1]);
    lifetime_norm=figure('units','normalized','outerposition',[0 0 1 1]);
%     lifetime_FZcorr=figure('units','normalized','outerposition',[0 0 1 1]);
%     lifetime_FZcorr_norm=figure('units','normalized','outerposition',[0 0 1 1]);
%     lifetime_mcSi_harm=figure('units','normalized','outerposition',[0 0 1 1]);
%     lifetime_mcSi_harm_norm=figure('units','normalized','outerposition',[0 0 1 1]);
%     lifetime_mcSi_ratio=figure('units','normalized','outerposition',[0 0 1 1]);
%     lifetime_mcSi_ratio_norm=figure('units','normalized','outerposition',[0 0 1 1]);
    Ntstar_time=figure('units','normalized','outerposition',[0 0 1 1]);
    labels_raw = {}; 
%     labels_FZcorr = {}; 
%     labels_mcSi = {};
%     hFZ = []; hFZnorm = [];
    cm = colormap(hsv(samp));
    for j = 1:samp
        index = find(strcmp(group_now{1,j},samples)==1);
        raw_now = lifetime_all{index}; 
        norm_now = norm_lifetime_all{index};
%         FZ_corr_now = FZ_corr{index}; 
%         FZ_corr_norm_now = FZ_corr_norm{index}; 
%         mcSi_harm_now = mcSi_harm{index}; 
%         mcSi_harm_norm_now = mcSi_harm_norm{index}; 
%         mcSi_ratio_now = mcSi_ratio{index}; 
%         mcSi_ratio_norm_now = mcSi_ratio_norm{index}; 
        Ntstar_now = Ntstar{index}; 
        %correct any 0 time starts
        if isempty(Ntstar_now)==0 && Ntstar_now(1,1) == 0
            Nstar_now(1,1) = 1; 
        end
        if isempty(raw_now)==0 && raw_now(1,1) == 0
            raw_now(1,1) = 1; 
        end
        if isempty(norm_now)==0 && norm_now(1,1) == 0 
            norm_now(1,1) = 1; 
        end
%         if isempty(FZ_corr_now)==0 && FZ_corr_now(1,1) == 0 
%             FZ_corr_now(1,1) = 1; 
%         end
%         if isempty(FZ_corr_norm_now)==0 && FZ_corr_norm_now(1,1) == 0 
%             FZ_corr_norm_now(1,1) = 1; 
%         end
%         if isempty(mcSi_harm_now)==0 && mcSi_harm_now(1,1) == 0 
%             mcSi_harm_now(1,1) = 1; 
%         end
%         if isempty(mcSi_harm_norm_now)==0 && mcSi_harm_norm_now(1,1) == 0 
%             mcSi_harm_norm_now(1,1) = 1; 
%         end
%         if isempty(mcSi_ratio_now)==0 && mcSi_ratio_now(1,1) == 0 
%             mcSi_ratio_now(1,1) = 1; 
%         end
%         if isempty(mcSi_ratio_norm_now)==0 && mcSi_ratio_norm_now(1,1) == 0 
%             mcSi_ratio_norm_now(1,1) = 1; 
%         end
        figure(Ntstar_time); 
        loglog(Ntstar_now(:,1),Ntstar_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
        hold all; 
        figure(lifetime_raw); 
        loglog(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
        hold all; 
        figure(lifetime_norm); 
        semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
        hold all;  
        labels_raw{end+1} = group_now{2,j};  
%         if strcmp(group_now{1,j},'FZ')==0 && strcmp(group_now{1,j},'FZ-new')==0 && ...
%             strcmp(group_now{1,j},'68-4')==0 && strcmp(group_now{1,j},'60a')==0 && ...
%             strcmp(group_now{1,j},'56b')==0 && strcmp(group_now{1,j},'FZ-new2')==0
%             for k = 1:length(index_FZ)
%                 figure(lifetime_FZcorr); 
%                 if k == 1
%                     hFZ(end+1)=loglog(FZ_corr_now(:,1),FZ_corr_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:)); 
%                     hold on; 
%                     figure(lifetime_FZcorr_norm); 
%                     hFZnorm(end+1)=semilogx(FZ_corr_norm_now(:,1),FZ_corr_norm_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:));
%                     hold on;  
%                 else
%                     loglog(FZ_corr_now(:,1),FZ_corr_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:)); 
%                     hold on; 
%                     figure(lifetime_FZcorr_norm); 
%                     semilogx(FZ_corr_norm_now(:,1),FZ_corr_norm_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:));
%                     hold on;  
%                 end
%             end
%             labels_FZcorr{end+1} = group_now{2,j}; 
%         end
%         if strcmp(group_now{1,j},'FZ')==0 && strcmp(group_now{1,j},'FZ-new')==0 && ...
%             strcmp(group_now{1,j},'FZ-12')==0 && strcmp(group_now{1,j},'66-2')==0 && ...
%             strcmp(group_now{1,j},'68-4')==0 && strcmp(group_now{1,j},'60a')==0 && ...
%             strcmp(group_now{1,j},'56b')==0 && strcmp(group_now{1,j},'C-1')==0 && ...
%             strcmp(group_now{1,j},'C-2')==0 && strcmp(group_now{1,j},'H-1')==0 && ...
%             strcmp(group_now{1,j},'H-2')==0 && strcmp(group_now{1,j},'FZ-new2')==0
%             figure(lifetime_mcSi_harm); 
%             loglog(mcSi_harm_now(:,1),mcSi_harm_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
%             hold all; 
%             figure(lifetime_mcSi_harm_norm); 
%             semilogx(mcSi_harm_norm_now(:,1),mcSi_harm_norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
%             hold all;
%             figure(lifetime_mcSi_ratio); 
%             semilogx(mcSi_ratio_now(:,1),mcSi_ratio_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
%             hold all; 
%             figure(lifetime_mcSi_ratio_norm); 
%             semilogx(mcSi_ratio_norm_now(:,1),mcSi_ratio_norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
%             hold all;
%             labels_mcSi{end+1} = group_now{2,j}; 
%         end
    end
    figure(lifetime_raw); 
    xlabel('time [s]','FontSize',25); 
    ylabel('lifetime [s]','FontSize',25); 
    legend(labels_raw); 
    xlim([1 max_time]); 
    title(plotting_names{i},'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_raw,[savedirname '\' plotting_names{i} '_raw' savename]);
    print(lifetime_raw,'-dpng','-r0',[savedirname '\' plotting_names{i} '_raw' savename '.png']);
    
    figure(Ntstar_time); 
    xlabel('time [s]','FontSize',25); 
    ylabel('N_t^* [s^-^1]','FontSize',25); 
    legend(labels_raw); 
    xlim([1 max_time]); 
    title(plotting_names{i},'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(Ntstar_time,[savedirname '\' plotting_names{i} '_Ntstar' savename]);
    print(Ntstar_time,'-dpng','-r0',[savedirname '\' plotting_names{i} '_Ntstar' savename '.png']);
    
    figure(lifetime_norm); 
    xlabel('time [s]','FontSize',25); 
    ylabel('normalized lifetime [-]','FontSize',25); 
    axis([1 max_time 0 2]);
    legend(labels_raw); 
    title(plotting_names{i},'FontSize',25); 
    set(0,'defaultAxesFontSize', 20)
    hgsave(lifetime_norm,[savedirname '\' plotting_names{i} '_norm' savename]);
    print(lifetime_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_norm' savename '.png']);
    
%     figure(lifetime_FZcorr); 
%     xlabel('time [s]','FontSize',25); 
%     ylabel('lifetime [s]','FontSize',25); 
%     legend(hFZ,labels_FZcorr); 
%     xlim([1 max_time]); 
%     title([plotting_names{i} ' corrected by FZ SRV'],'FontSize',25); 
%     set(0,'defaultAxesFontSize', 20)
%     hgsave(lifetime_FZcorr,[savedirname '\' plotting_names{i} '_FZcorr' savename]);
%     print(lifetime_FZcorr,'-dpng','-r0',[savedirname '\' plotting_names{i} '_FZcorr' savename '.png']);
%     
%     figure(lifetime_FZcorr_norm); 
%     xlabel('time [s]','FontSize',25); 
%     ylabel('normalized lifetime [-]','FontSize',25); 
%     axis([1 max_time 0 2]);
%     legend(hFZnorm,labels_FZcorr); 
%     title([plotting_names{i} ' corrected by FZ SRV'],'FontSize',25); 
%     set(0,'defaultAxesFontSize', 20)
%     hgsave(lifetime_FZcorr_norm,[savedirname '\' plotting_names{i} '_FZcorr_norm' savename]);
%     print(lifetime_FZcorr_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_FZcorr_norm' savename '.png']);
%     
%     figure(lifetime_mcSi_harm); 
%     xlabel('time [s]','FontSize',25); 
%     ylabel('lifetime [s]','FontSize',25); 
%     legend(labels_mcSi); 
%     xlim([1 max_time]); 
%     title([plotting_names{i} ' corrected by mcSi control harmonic sum'],'FontSize',25); 
%     set(0,'defaultAxesFontSize', 20)
%     hgsave(lifetime_mcSi_harm,[savedirname '\' plotting_names{i} '_mcSi_harm' savename]);
%     print(lifetime_mcSi_harm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_harm' savename '.png']);
%     
%     figure(lifetime_mcSi_harm_norm); 
%     xlabel('time [s]','FontSize',25); 
%     ylabel('normalized lifetime [-]','FontSize',25); 
%     axis([1 max_time 0 2]);
%     legend(labels_mcSi); 
%     title([plotting_names{i} ' corrected by mcSi control harmonic sum'],'FontSize',25); 
%     set(0,'defaultAxesFontSize', 20)
%     hgsave(lifetime_mcSi_harm_norm,[savedirname '\' plotting_names{i} '_mcSi_harm_norm' savename]);
%     print(lifetime_mcSi_harm_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_harm_norm' savename '.png']);
%     
%     figure(lifetime_mcSi_ratio); 
%     xlabel('time [s]','FontSize',25); 
%     ylabel('lifetime [s]','FontSize',25); 
%     legend(labels_mcSi); 
%     xlim([1 max_time]);
%     title([plotting_names{i} ' corrected by mcSi control ratio'],'FontSize',25); 
%     set(0,'defaultAxesFontSize', 20)
%     hgsave(lifetime_mcSi_ratio,[savedirname '\' plotting_names{i} '_mcSi_ratio' savename]);
%     print(lifetime_mcSi_ratio,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_ratio' savename '.png']);
%     
%     figure(lifetime_mcSi_ratio_norm); 
%     xlabel('time [s]','FontSize',25); 
%     ylabel('normalized lifetime [-]','FontSize',25); 
%     axis([1 max_time 0 2]);
%     legend(labels_mcSi); 
%     title([plotting_names{i} ' corrected by mcSi control ratio'],'FontSize',25); 
%     set(0,'defaultAxesFontSize', 20)
%     hgsave(lifetime_mcSi_ratio_norm,[savedirname '\' plotting_names{i} '_mcSi_ratio_norm' savename]);
%     print(lifetime_mcSi_ratio_norm,'-dpng','-r0',[savedirname '\' plotting_names{i} '_mcSi_ratio_norm' savename '.png']);
end

%Let's save all of the processed data as a .mat file for easy use later
% save([savedirname '\' bora '_all_data' savename '.mat'],'lifetime_all',...
%     'norm_lifetime_all','FZ_corr','FZ_corr_norm','mcSi_harm',...
%     'mcSi_harm_norm','mcSi_ratio','mcSi_ratio_norm','doping_all',...
%     'thickness_all','samples','SRV_FZ','surface_control','Ntstar');
save([savedirname '\' bora '_all_data' savename '.mat'],'lifetime_all',...
    'norm_lifetime_all','doping_all',...
    'thickness_all','samples','Ntstar');

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
        %correct any 0 time starts
        if isempty(raw_now)==0 && raw_now(1,1) == 0
            raw_now(1,1) = 1; 
        end
        if isempty(norm_now)==0 && norm_now(1,1) == 0 
            norm_now(1,1) = 1; 
        end
        if isempty(FZ_corr_now)==0 && FZ_corr_now(1,1) == 0 
            FZ_corr_now(1,1) = 1; 
        end
        if isempty(FZ_corr_norm_now)==0 && FZ_corr_norm_now(1,1) == 0 
            FZ_corr_norm_now(1,1) = 1; 
        end
        if isempty(mcSi_harm_now)==0 && mcSi_harm_now(1,1) == 0 
            mcSi_harm_now(1,1) = 1; 
        end
        if isempty(mcSi_harm_norm_now)==0 && mcSi_harm_norm_now(1,1) == 0 
            mcSi_harm_norm_now(1,1) = 1; 
        end
        if isempty(mcSi_ratio_now)==0 && mcSi_ratio_now(1,1) == 0 
            mcSi_ratio_now(1,1) = 1; 
        end
        if isempty(mcSi_ratio_norm_now)==0 && mcSi_ratio_norm_now(1,1) == 0 
            mcSi_ratio_norm_now(1,1) = 1; 
        end
        figure(lifetime_raw); 
        loglog(raw_now(:,1),raw_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
        hold all; 
        figure(lifetime_norm); 
        semilogx(norm_now(:,1),norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
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
                    hFZ(end+1)=loglog(FZ_corr_now(:,1),FZ_corr_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:)); 
                    hold on; 
                    figure(lifetime_FZcorr_norm); 
                    hFZnorm(end+1)=semilogx(FZ_corr_norm_now(:,1),FZ_corr_norm_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:));
                    hold on;  
                else
                    loglog(FZ_corr_now(:,1),FZ_corr_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:)); 
                    hold on; 
                    figure(lifetime_FZcorr_norm); 
                    semilogx(FZ_corr_norm_now(:,1),FZ_corr_norm_now(:,k+1),'-o','LineWidth',3,'MarkerSize',10,'color',cm(j,:));
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
            loglog(mcSi_harm_now(:,1),mcSi_harm_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
            hold all; 
            figure(lifetime_mcSi_harm_norm); 
            semilogx(mcSi_harm_norm_now(:,1),mcSi_harm_norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
            hold all;
            figure(lifetime_mcSi_ratio); 
            semilogx(mcSi_ratio_now(:,1),mcSi_ratio_now(:,2),'-o','LineWidth',3,'MarkerSize',10); 
            hold all; 
            figure(lifetime_mcSi_ratio_norm); 
            semilogx(mcSi_ratio_norm_now(:,1),mcSi_ratio_norm_now(:,2),'-o','LineWidth',3,'MarkerSize',10);
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
close all; clc; 
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

%loop over all the samples and get the relevant data in a structure we can
%work with
for i = 1:length(samples)
    %figure out which measurements were taking for this sample
    meas_thissample = meas(i,:); 
    lifetime_store = {};
    doping_store = [];
    thickness_store = [];
    %We will do the same procedure for each measurement at each degradation
    %time
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
lifetime_breakdown = figure; 
defect_summary = cell(num_samples,1); 
for i = 1:num_samples
    index = find(strcmp(samples,lifetime_analysis{1,i})==1);
    raw_now = lifetime_all{index,2}; 
    raw_times = lifetime_all{index,1}; 
    [num_meas,columns] = size(raw_now); 
    doping_now = doping_all{index};
    thickness_now = thickness_all{index}; 
    %Check that the parameters are consistent for this sample
    if length(doping_now) ~= length(find(doping_now(1)==doping_now))
        disp(['the doping for sample ' samples{index} ' is not consistent']); 
    end
    if length(thickness_now) ~= length(find(thickness_now(1)==thickness_now))
        disp(['the thickness for sample ' samples{index} ' is not consistent']); 
    end
    if strcmp(lifetime_analysis{1,i},'FZ')==0 && strcmp(lifetime_analysis{1,i},'FZ-new')==0 && ...
            strcmp(lifetime_analysis{1,i},'68-4')==0 && strcmp(lifetime_analysis{1,i},'60a')==0 && ...
            strcmp(lifetime_analysis{1,i},'56b')==0 && strcmp(lifetime_analysis{1,i},'FZ-new2')==0
        %Then we explicitly find the lifetime using the FZ wafer as a
        %reference
        tau_rev = zeros(num_meas,length(index_FZ)); 
        %We'll want to save the fit data in a structure that makes sense
        defect_fits = cell(num_meas,length(index_FZ)); 
        for k = 1:length(index_FZ)
            SRV_now = SRV_FZ{k,2};
            SRV_times = SRV_FZ{k,1}; 
            count = 1; 
            for j = 1:num_meas
                %Get the length of the measured lifetime
                lifetime_now = raw_now{j}; 
                [injections,columns] = size(lifetime_now); 
                %Figure out which SRV matters here
                try
                    t_index = find(SRV_times(:,1)==raw_times(j,1));
                    SRV_this_time = SRV_now{t_index}; 
                    %interpolate the SRV at the measured injection levels
                    [SRV] = interp1(SRV_this_time(:,1),SRV_this_time(:,2),lifetime_now(:,1)); 
                    %Get the doping and thickness of this sample
                    W_now = thickness_now(j); 
                    Na_now = doping_now(j); 
                    %loop over the injection levels and calculate diffusivity,
                    %surface lifetime, intrinsic lifetime to get at the
                    %injection-dependent SRH lifetime. 
                    tau_surf = zeros(injections,1); 
                    tau_intr = zeros(injections,1); 
                    for m = 1:injections
                        [D_now,Dh] = diffusivity(300,'p',Na_now,lifetime_now(m,1));
                        tau_surf(m) = (W_now./(2.*SRV(m)))+((1/D_now).*((W_now/pi)^2)); %cm/s
                        tau_intr(m) = Richter(300,lifetime_now(m,1),Na_now,'p');
                    end
                    %Let's simplify the surface lifetime - it should be
                    %basically constant
                    surf_cutoff_high = 2e15; 
                    [deltan_surf_rev,tau_surf_rev]=...
                        remove_highinj(lifetime_now(:,1),tau_surf,surf_cutoff_high);
                    surf_avg = nanmean(tau_surf_rev); 
                    tau_surf_new = ones(size(tau_intr)).*surf_avg; 
                    tau_SRH = ((1./lifetime_now(:,2))-(1./tau_surf_new)-(1./tau_intr)).^(-1);
                    figure(lifetime_breakdown); clf; 
                    loglog(lifetime_now(:,1),lifetime_now(:,2)); 
                    hold all; 
                    loglog(lifetime_now(:,1),tau_surf); 
                    hold all;
                    loglog(lifetime_now(:,1),tau_surf_new); 
                    hold all;
                    loglog(lifetime_now(:,1),tau_intr); 
                    hold all;
                    loglog(lifetime_now(:,1),tau_SRH); 
                    legend('measured','surface','simple surface','intrinsic','SRH'); 
                    if count == 1
                        disp('Select the region for cutting off the HIGH injection data');
                        [cutoff_high,nothing]=ginput(1);
                        disp('Select the region for cutting off the LOW injection data');
                        [cutoff_low,nothing]=ginput(1);
                        count = 2; 
                    end
                    mkdir([savedirname '\' lifetime_analysis{1,i}]);
                    save_this = [savedirname '\' lifetime_analysis{1,i} '\' ...
                        lifetime_analysis{1,i} '_' num2str(raw_times(j,1)) ...
                        's_' surface_control{k}];
                    %get rid of any NaN values first
                    index_nan = find(isnan(tau_SRH)==1); 
                    tau_SRH(index_nan)=[]; 
                    deltan_SRH = lifetime_now(:,1); deltan_SRH(index_nan)=[];
                    [easy_summary,all_defect] = fit_procedure(lifetime_breakdown,...
                        deltan_SRH,tau_SRH,save_this,300,Na_now,type,...
                        cutoff_low,cutoff_high);
                    easy_summary = [raw_times(j,1) easy_summary]; 
                    defect_fits{j,k} = {easy_summary,all_defect}; 
                catch
                    %This time does not exist for this FZ wafer
                    defect_fits{j,k} = []; 
                end
            end
        end
    elseif strcmp(lifetime_analysis{1,i},'68-4')==1 || ...
            strcmp(lifetime_analysis{1,i},'60a')==1 || ...
            strcmp(lifetime_analysis{1,i},'56b')==1
        %Then we will go ahead and use the harmonic sum 
        %Set the initial data
        raw_initial = raw_now{1};
        deltan_initial = raw_initial(:,1); 
        tau_initial = raw_initial(:,2); 
        defect_fits = cell(num_meas-1,1); 
        count = 1; 
        %Loop over all measured data
        for j = 2:num_meas
            W_now = thickness_now(j); 
            Na_now = doping_now(j);
            raw_this_time = raw_now{j}; 
            deltan = raw_this_time(:,1); 
            tau = raw_this_time(:,2); 
            figure(lifetime_breakdown); clf; 
            loglog(deltan_initial,tau_initial); 
            hold all; 
            loglog(deltan,tau); 
            hold all; 
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
            tau_SRH = ((1./tau)-(1./tau_initial)).^(-1);
            loglog(deltan,tau_SRH);
            legend('initial','this time','SRH'); 
            if count == 1
                disp('Select the region for cutting off the HIGH injection data');
                [cutoff_high,nothing]=ginput(1);
                disp('Select the region for cutting off the LOW injection data');
                [cutoff_low,nothing]=ginput(1);
                count = 2; 
            end
            mkdir([savedirname '\' lifetime_analysis{1,i}]);
            save_this = [savedirname '\' lifetime_analysis{1,i} '\' ...
                lifetime_analysis{1,i} '_' num2str(raw_times(j,1)) ...
                's_harmSum'];
            [easy_summary,all_defect] = fit_procedure(lifetime_breakdown,...
                deltan,tau_SRH,save_this,300,Na_now,type,cutoff_low,cutoff_high);
            easy_summary = [raw_times(j,1) easy_summary]; 
            defect_fits{j,k} = {easy_summary,all_defect}; 
        end
    end
    defect_summary{i} = defect_fits; 
end

save([savedirname '\defect_fits.mat'],'defect_summary','lifetime_analysis','surface_control'); 
%% Look through the data from the previous fit
close all; clc; 
load([savedirname '\defect_fits.mat']); 
load([dirname '\' bora '_all_data' savename '.mat']); 
%Loop over all of the investigated samples
[rows,num_samples] = size(lifetime_analysis); 
for i = 1:num_samples
    defect_now = defect_summary{i}; 
    [num_meas,columns] = size(defect_now); 
    time_vector = zeros(num_meas,1); 
    k1 = zeros(num_meas,1); 
    k2 = zeros(num_meas,1); 
    count = 1; 
    %We will plot the E-k curves together for defects 1,2 for better
    %visibility
    allEk=figure('units','normalized','outerposition',[0 0 1 1]);
    labels = {};
    for j = 1:num_meas
        %Figure out which FZ control samples mattered
        indexFZ = find(~cellfun(@isempty,defect_now(j,:)));
        for k = 1:length(indexFZ)
            all_defect = defect_now{j,indexFZ(k)}; 
            easy_now = all_defect{1}; 
            if easy_now(1) == 0
                time_vector(count) = 1; 
            else
                time_vector(count) = easy_now(1); 
            end
            k1(count) = easy_now(2); 
            k2(count) = easy_now(4); 
                        
            %The above was the easy summary
            %Now we do the more challenging summary - details for each
            %sample
            full_details = all_defect{2}; 
            figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(3,2,[1 2]);
            %linearized lifetime
            semilogy(full_details{6},full_details{7},'k','LineWidth',3); 
            def1 = (full_details{1}(1,1).*full_details{6})+full_details{1}(1,2);
            def2 = (full_details{1}(2,1).*full_details{6})+full_details{1}(2,2);
            hold all;
            semilogy(full_details{6},def1,'r--','LineWidth',3); 
            hold all;
            semilogy(full_details{6},def2,'b--','LineWidth',3); 
            together = ((1./def1)+(1./def2)).^(-1); 
            hold all;
            plot(full_details{6},together,'x'); 
            xlabel('X = n/p [-]'); 
            ylabel('\tau_{SRH} [s]'); 
            title([lifetime_analysis{1,i} ', ' lifetime_analysis{2,i} ...
                ', ' num2str(time_vector(count)) 's']); 
            %defect 1 k value
            subplot(3,2,3); 
            semilogy(full_details{3}{1},full_details{4}{1},'r','LineWidth',3); 
            xlabel('E_t-E_i [eV]'); 
            ylabel('k [-]'); 
            xlim([min(full_details{3}{1}) max(full_details{3}{1})]); 
            title('default defect 1'); 
            %defect 1 tau_n0
            subplot(3,2,5); 
            semilogy(full_details{3}{1},1./full_details{5}{1},'r','LineWidth',3); 
            xlabel('E_t-E_i [eV]'); 
            ylabel('\tau_{n0} [s]'); 
            xlim([min(full_details{3}{1}) max(full_details{3}{1})]);
            %defect 2 k value
            subplot(3,2,4); 
            semilogy(full_details{3}{2},full_details{4}{2},'b','LineWidth',3); 
            xlabel('E_t-E_i [eV]'); 
            ylabel('k [-]'); 
            xlim([min(full_details{3}{1}) max(full_details{3}{1})]);
            title('default defect 2'); 
            %defect 2 tau_n0
            subplot(3,2,6);
            semilogy(full_details{3}{2},1./full_details{5}{2},'b','LineWidth',3); 
            xlabel('E_t-E_i [eV]'); 
            ylabel('\tau_{n0} [s]'); 
            xlim([min(full_details{3}{1}) max(full_details{3}{1})]);
            save_this = [savedirname '\' lifetime_analysis{1,i} '\' ...
                lifetime_analysis{1,i} '_' num2str(count) '_defect_details'];
            hgsave(gcf,save_this);
            print(gcf,'-dpng','-r0',[save_this '.png']);
            
            figure(allEk); 
            subplot(2,2,1); 
            semilogy(full_details{3}{1},full_details{4}{1},'LineWidth',3);
            hold all; 
            subplot(2,2,2);
            semilogy(full_details{3}{2},full_details{4}{2},'LineWidth',3);
            hold all; 
            subplot(2,2,3); 
            semilogy(full_details{3}{1},1./full_details{5}{1},'LineWidth',3); 
            hold all; 
            subplot(2,2,4);
            semilogy(full_details{3}{2},1./full_details{5}{2},'LineWidth',3); 
            hold all; 
            labels{count} = [num2str(time_vector(count)) 's'];
            count = count+1; 
        end
    end
    %find the lifetime degradation data
    index_sample = find(strcmp(samples,lifetime_analysis{1,i})==1); 
    norm_now = norm_lifetime_all{index_sample}; 
    if norm_now(1,1) == 0
        norm_now(1,1) = 1; 
    end
    %one figure per sample - summary of all times
    figure('units','normalized','outerposition',[0 0 1 1]);
    %top plot - lifetime degradation
    subplot(3,1,1); 
    semilogx(norm_now(:,1),norm_now(:,2),'LineWidth',3); 
    ylabel('normalized lifetime [-]','FontSize',14); 
    axis([1 max_time*1.25 0 1.5]);
    title(lifetime_analysis{2,i},'FontSize',20); 
    %middle plot - k1 defect at midgap
    subplot(3,1,2); 
    semilogx(time_vector,k1,'-','LineWidth',3); 
    ylabel('k [-]','FontSize',14); 
    xlim([1 max_time*1.25]);
    title('dominant defect','FontSize',14); 
    %bottom plot - k2 defect at midgap
    subplot(3,1,3); 
    semilogx(time_vector,k2,'-','LineWidth',3); 
    ylabel('k [-]','FontSize',14); 
    xlabel('degradation time [s]','FontSize',14); 
    xlim([1 max_time*1.25]);
    title('secondary defect','FontSize',14); 
    save_this = [savedirname '\' lifetime_analysis{1,i} '\' ...
                lifetime_analysis{1,i} '_defect_summary'];
    hgsave(gcf,save_this);
    print(gcf,'-dpng','-r0',[save_this '.png']);
    
    figure(allEk); 
    subplot(2,2,1);
    xlabel('E_t-E_i [eV]'); 
    ylabel('k [-]'); 
    xlim([min(full_details{3}{1}) max(full_details{3}{1})]); 
    legend(labels); 
    title('default defect 1'); 
    subplot(2,2,2); 
    xlabel('E_t-E_i [eV]'); 
    ylabel('k [-]'); 
    xlim([min(full_details{3}{1}) max(full_details{3}{1})]);
    title('default defect 2');
    subplot(2,2,3); 
    xlabel('E_t-E_i [eV]'); 
    ylabel('\tau_{n0} [s]'); 
    xlim([min(full_details{3}{1}) max(full_details{3}{1})]);
    subplot(2,2,4); 
    xlabel('E_t-E_i [eV]'); 
    ylabel('\tau_{n0} [s]'); 
    xlim([min(full_details{3}{1}) max(full_details{3}{1})]);
    save_this = [savedirname '\' lifetime_analysis{1,i} '\' ...
                lifetime_analysis{1,i} '_all_E_k_taun0'];
    hgsave(gcf,save_this);
    print(gcf,'-dpng','-r0',[save_this '.png']);
               
    close all; 
end

%% Look through data already calculated and extract data for further analysis
close all; clc; 
load([savedirname '\defect_fits.mat']); 
load([dirname '\' bora '_all_data' savename '.mat']); 
%Loop over all of the investigated samples
[rows,num_samples] = size(lifetime_analysis); 
for i = 1:num_samples
    %File for saving the full injection dependent data for fitting
    tau_file = [savedirname '\' lifetime_analysis{1,i} '\IDLS two and three curve fitting.xlsm']; 
    defect_now = defect_summary{i}; 
    [num_meas,columns] = size(defect_now); 
    %Index for referencing the original data
    index = find(strcmp(samples,lifetime_analysis{1,i})==1);
    raw_tau = lifetime_all{index}; 
    norm_raw = norm_lifetime_all{index}; 
    Ntstar_raw = Ntstar{index}; 
    %We want to write 7 columns of summary data for plotting
    to_write = zeros(num_meas,7); 
    FZ_track = cell(num_meas,1); 
    count = 1; 
    if strcmp(lifetime_analysis{1,i},'68-4')==1 || ...
        strcmp(lifetime_analysis{1,i},'60a')==1 || ...
        strcmp(lifetime_analysis{1,i},'56b')==1
        %We need to make up for the first measurement
        if raw_tau(1,1) == 0
        	to_write(count,1) = 1;
        else
            to_write(count,1) = raw_tau(1,1); 
        end
        to_write(count,2) = raw_tau(1,2);
        to_write(count,4) = norm_raw(1,2); 
        to_write(count,6) = Ntstar_raw(1,2);
        count = count+1; 
    end
    for j = 1:num_meas
        %Figure out which FZ control samples mattered
        indexFZ = find(~cellfun(@isempty,defect_now(j,:)));
        for k = 1:length(indexFZ)
            all_defect = defect_now{j,indexFZ(k)}; 
            easy_now = all_defect{1}; 
            %write the time
            if easy_now(1) == 0
                to_write(count,1) = 1; 
            else
                to_write(count,1) = easy_now(1); 
            end
            %write the raw lifetime, we'll assume for now that what we
            %previously processed is at the right injection level
            index_time = find(raw_tau(:,1)==easy_now(1));
            if isempty(index_time)==0
                to_write(count,2) = raw_tau(index_time,2); 
            end
            %write the SRH lifetime
            tau_SRH = all_defect{2}{7}; 
            deltan_SRH = all_defect{2}{8}; 
            X_SRH = all_defect{2}{6}; 
            tau_SRH_single = interp1(deltan_SRH,tau_SRH,deltan_target); 
            to_write(count,3) = tau_SRH_single; 
            %write the normalized raw lifetime
            index_time = find(norm_raw(:,1)==easy_now(1));
            if isempty(index_time)==0
                to_write(count,4) = norm_raw(index_time,2); 
            end
            %write the Ntstar from the raw lifetime
            index_time = find(Ntstar_raw(:,1)==easy_now(1));
            if isempty(index_time)==0
                to_write(count,6) = Ntstar_raw(index_time,2); 
            end
            
            if strcmp(lifetime_analysis{1,i},'68-4')==1 || ...
                strcmp(lifetime_analysis{1,i},'60a')==1 || ...
                strcmp(lifetime_analysis{1,i},'56b')==1
                FZ_track{count} = 'harmonic'; 
                if count == 2
                    to_write(count,5) = 1; 
                else
                    to_write(count,5) = tau_SRH_single/to_write(2,2); 
                end
                %write the Ntstar from the SRH lifetime
                to_write(count,7) = (1/to_write(count,3))-(1/to_write(2,3)); 
            else
                FZ_track{count} = surface_control{indexFZ};
                %write the normalized SRH lifetime
                if count == 1
                    to_write(count,5) = 1; 
                else
                    to_write(count,5) = tau_SRH_single/to_write(1,2);
                end
                %write the Ntstar from the SRH lifetime
                to_write(count,7) = (1/to_write(count,3))-(1/to_write(1,3)); 
            end
            %write the injection dependent data for fit refining
            tau_write = [X_SRH,tau_SRH]; 
            Xs = length(X_SRH); 
            xlsrange = ['B9:C' num2str(9+Xs)]; 
            xlswrite(tau_file,tau_write,num2str(count),xlsrange); 
            %also write the starting fit parameters for two defects
            p_write = all_defect{2}{1}; 
            p_write = [p_write(1,1);p_write(1,2);p_write(2,1);p_write(2,2)];
            xlswrite(tau_file,p_write,num2str(count),'L2:L5');
            count = count+1; 
        end
    end
    save_this = [savedirname '\' lifetime_analysis{1,i} '\all_data_summary.xlsx'];
    labels = {'time [s]','tau [s]','tauSRH [s]','norm tau [-]','norm tauSRH [-]','Nt*','Nt* SRH','surface'}; 
    xlswrite(save_this,labels,1,'A1:H1'); 
    xlswrite(save_this,to_write,1,'A2'); 
    xlswrite(save_this,FZ_track,1,'H2'); 
end   
    
