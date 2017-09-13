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

%This script analyzes lifetime measurements taken with WCT-120. 
clear all; close all; 

%Location of your data
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\PERC LeTID Advanced System\Files after measurement\50C'; 
%Where you want the data to save (could be same as directory)
processing_directory = 'C:\Users\Mallory Jensen\Documents\LeTID\PERC LeTID Advanced System\Files after measurement\50C';
%Where the raw data for the float zone control is
SRV_directory = 'C:\Users\Malloryj\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\By sample\16-6-28-P-2\Summary files';
%type - p or n
type = 'n'; 
%% Collect and process raw data - WCT-120 files
%Find all of the files in the directory
[fileList,fileListShort] = getAllFiles(directory); 
%Get the injection-dependent lifetime data for all files
process_xls_data(directory,[directory '\Raw_data.mat']);
%We also want to store all of the information for each file
%T, thickness, resistivity (entered/measured), type, optical constant, calibration,
%1/64 or 1/1
for file = 1:length(fileList)
    this_file = fileList{file};
    thick{file,1} = xlsread(this_file,'User','B6');
    res{file,1} = xlsread(this_file,'User','C6');
    oc{file,1} = xlsread(this_file,'User','E6');
    temp{file,1} = 25;
    meas_res{file,1} = xlsread(this_file,'Summary','N2');
    calib{file,1} = xlsread(this_file,'Settings','C5');
    doping{file,1} = xlsread(this_file,'Summary','E2');
end
info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
save([directory '\meas_info.mat'],'info'); 

%% Now that the raw data has been processed, plot the curves together
load([directory '\Raw_data.mat']);
load([directory '\meas_info.mat']); 
h=figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(dataSave)
    dataNow = dataSave{i}; 
    curves(i)=loglog(dataNow(:,1),dataNow(:,2),'.'); 
    hold all; 
    label{i} = [info(i).filename];
end
%% Label and save the current summary plot
xlabel('Excess carrier density [cm^-^3]','FontSize',20); 
ylabel('Lifetime [s]','FontSize',20);
legend(curves',label');
hgsave(h,[directory '\Lifetime w T']);
print(h,'-dpng','-r0',[directory '\Lifetime w T.png']); 
%% Load the data that will be processed and lifetime contributions - no surface contribution
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    for j = 1:length(deltan_rev)
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
    end
    %Finally, we can calculate the SRH lifetime, assuming no surface
    %contribution
    tau_SRH = ((1./tau_rev)-(1./tau_intr)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6);
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6);
    hold all;
    loglog(deltan_rev,tau_SRH.*1e6);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Intrinsic','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');

%% Load the data and calculate the SRV as a function of injection for each curve
%First load the data
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']);
SRVtoSave = cell(length(dataSave),1);
h=figure; 
label = zeros(length(dataSave),1)
%Loop through each temperature calculate the SRV
for i = 1:length(dataSave)
    dataNow = dataSave{i}; 
    tau = dataNow(:,2); 
    deltan = dataNow(:,1); 
    %We might want to trim the lifetime
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %We need the intrinsic lifetime
    deltan = deltan_rev;
    tau = tau_rev; 
    tau_intr = zeros(size(deltan));
    diffusivity_save = zeros(size(deltan)); 
    for j = 1:length(deltan)
        %Get the intrinsic lifetime
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan(j),info(i).doping,type);
        %unfortunately the diffusivity is also temperature and injection
        %dependent
        [De,Dh] = diffusivity(info(i).temperature+273.15,type,info(i).doping,deltan(j));
        if type == 'n'
            diffusivity_save(j,1) = Dh; %hole is minority carrier
        elseif type == 'p'
            diffusivity_save(j,1) = De; %electron is minority carrier
        end
    end
    %Calculate the surface-related lifetime assuming zero SRH contribution
    tau_surf = ((1./tau)-(1./tau_intr)).^(-1);
    %Calculate the SRV, including the injection-dependent diffusivity
    SRV = info(i).thickness./((tau_surf-((1./diffusivity_save).*((info(i).thickness/pi)^2))).*2);
    SRVtoSave{i} = [deltan SRV]; 
    %Plot the result so that we can compare
    figure(h);
    loglog(deltan,SRV,'LineWidth',2);
    label(i,1)=info(i).temperature; 
    hold all; 
end
xlabel('Excess carrier density [cm^-^3]','FontSize',20); 
ylabel('Surface recombination velocity [cm/s]','FontSize',20);
legend(num2str(label));
hgsave(h,[directory '\SRV w T']);
print(h,'-dpng','-r0',[directory '\SRV w T.png']); 
save([directory '\SRV_data.mat'],'SRVtoSave'); 

%% Load the data and process the lifetime contributions with SRV given by a .mat file
%Load the data for processing
load([processing_directory '\Raw_data.mat']);
load([processing_directory '\meas_info.mat']); 
%Load the SRV data
load([SRV_directory '\SRV_data.mat']); 
%Now, for each temperature we will be doing the same operations. Loop
%through. 
for i = 1:length(dataSave);
    dataNow = dataSave{i}; 
    tau = dataNow(:,2);
    deltan = dataNow(:,1);
    %We first need to remove the high injection data
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
    title(['Temperature = ' num2str(info(i).temperature)]);
    %We will also need the intrinsic lifetime
    tau_intr = zeros(size(deltan_rev));
    diffusivity_save = zeros(size(deltan)); 
    for j = 1:length(deltan_rev)
        %Get the intrinsic lifetime
        tau_intr(j,1) = Richter(info(i).temperature+273.15,deltan_rev(j),info(i).doping,type);
        %unfortunately the diffusivity is also temperature and injection
        %dependent
        [De,Dh] = diffusivity(info(i).temperature+273.15,type,info(i).doping,deltan_rev(j));
        if type == 'n'
            diffusivity_save(j,1) = Dh; %hole is minority carrier
        elseif type == 'p'
            diffusivity_save(j,1) = De; %electron is minority carrier
        end
    end
    %We need the surface contribution which comes from a reference wafer
    SRV_now = SRVtoSave{i};
    %Interpolate the SRV so that it matches the measured lifetime
    SRVq = interp1(SRV_now(:,1),SRV_now(:,2),deltan_rev); 
    tau_surf =(info(i).thickness./(2.*SRVq))+((1./diffusivity_save(i)).*((info(i).thickness/pi)^2)); %cm/s
    %Finally, we can calculate the SRH lifetime
    tau_SRH = ((1./tau_rev)-(1./tau_intr)-(1./tau_surf)).^(-1);
    %We should plot all of the contributions together
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6,'LineWidth',2);
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6,'LineWidth',2);
    hold all;
    loglog(deltan_rev,tau_surf.*1e6,'LineWidth',2); 
    hold all; 
    loglog(deltan_rev,tau_SRH.*1e6,'LineWidth',2);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Intrinsic','Surface','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    hgsave(h,[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature))]);
    print(h,'-dpng','-r0',[processing_directory '\Lifetime breakdown ' num2str(round(info(i).temperature)) '.png']); 
    %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev;
    tau_SRH_store{i,1} = tau_SRH;
    tau_intr_store{i,1} = tau_intr; 
    tau_surf_store{i,1} = tau_surf; 
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'tau_surf',tau_surf_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');

%% Load the data and ask for where to crop the SRH data based on the contributions
%Load the data for processing
load([processing_directory '\lifetime_breakdown.mat']);
for i = 1:length(dataSave)
    deltan_rev = lifetime_breakdown(i).deltan; 
    tau_rev = lifetime_breakdown(i).tau; 
    tau_SRH = lifetime_breakdown(i).tau_SRH; 
    tau_intr = lifetime_breakdown(i).tau_intr;
    tau_surf = lifetime_breakdown(i).tau_surf;
    figure;
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_rev.*1e6,'LineWidth',2);
    hold all; 
    loglog(deltan_rev,tau_intr.*1e6,'LineWidth',2);
    hold all;
    loglog(deltan_rev,tau_surf.*1e6,'LineWidth',2); 
    hold all; 
    loglog(deltan_rev,tau_SRH.*1e6,'LineWidth',2);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    legend('Measured','Intrinsic','Surface','SRH');
    title(['Temperature = ' num2str(info(i).temperature)]);
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_highinj(deltan_rev,tau_SRH,cutoff);
    [deltan_rev,tau_rev_rev] = remove_highinj(deltan_rev,tau_rev,cutoff);
    [deltan_rev,tau_surf_rev] = remove_highinj(deltan_rev,tau_surf,cutoff);
    [deltan_rev,tau_intr_rev] = remove_highinj(deltan_rev,tau_intr,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_lowinj(deltan_rev,tau_SRH_rev,cutoff);
    [deltan_rev,tau_rev_rev] = remove_lowinj(deltan_rev,tau_rev_rev,cutoff);
    [deltan_rev,tau_surf_rev] = remove_lowinj(deltan_rev,tau_surf_rev,cutoff);
    [deltan_rev,tau_intr_rev] = remove_lowinj(deltan_rev,tau_intr_rev,cutoff);
     %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_store{i,1} = tau_rev_rev;
    tau_SRH_store{i,1} = tau_SRH_rev;
    tau_intr_store{i,1} = tau_intr_rev; 
    tau_surf_store{i,1} = tau_surf_rev;
end
lifetime_breakdown = struct('tau',tau_store,'deltan',deltan_store,'tau_SRH',tau_SRH_store,'tau_intr',tau_intr_store,'tau_surf',tau_surf_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
%% Perform the lifetime fitting
%Load the data which has been separated into lifeitme contributions
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']); 
fit_tries = 1e5; 
%This loop may be slow due to fitting
for i = 1:length(lifetime_breakdown)
    %Given the SRH lifetime, we just need to linearize the injection level
    %Plot the Murphy linearization for n-type
    %Get sample parameters at specified temperature
    [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(info(i).temperature+273.15,info(i).doping,type); 
    %Normalized carrier density
    if type == 'p'
        X = (n0+lifetime_breakdown(i).deltan)./(p0+lifetime_breakdown(i).deltan);
    elseif type == 'n'
        X = (p0+lifetime_breakdown(i).deltan)./(n0+lifetime_breakdown(i).deltan);
    end
    tau_SRH_now = lifetime_breakdown(i).tau_SRH; 
    indices = find(tau_SRH_now<0 | isnan(tau_SRH_now)==1); 
    tau_SRH_now(indices) = [];
    X(indices) = [];
    X_store{i} = X; 
    tau_SRH_store{i} = tau_SRH_now;
    %Write the data to an XLS file for easy copying
    
    xlswrite([processing_directory '\Linearized_data.xlsx'],[X,tau_SRH_now],['Sheet' num2str(i)]); 
    %Now we fit. Fit in microseconds since I'm familiar with good MSE for
    %those units
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_SRH_now.*1e6,info(i).temperature,processing_directory,fit_tries);
    to_write = zeros(6,3); 
    to_write(1:2,1) = one_defect{i,1}';
    twodef = two_defects{i,1}; 
    to_write(1:2,2) = twodef(1,:)';
    to_write(3:4,2) = twodef(2,:)';
    threedef = three_defects{i,1}; 
    to_write(1:2,3) = threedef(1,:)';
    to_write(3:4,3) =threedef(2,:)';
    to_write(5:6,3) = threedef(3,:)';
    xlswrite([processing_directory '\Linearized_data.xlsx'],to_write,['Sheet' num2str(i)],'C1:E6'); 
end
best_fits = struct('one_defect',one_defect,'MSE_one',MSE_one,'two_defects',two_defects,'MSE_two',MSE_two,'three_defects',three_defects,'MSE_three',MSE_three,'all_fits',all_parameters_store,'all_MSE',all_MSE_store);
save([processing_directory '\best_fits.mat'],'best_fits');
save([processing_directory '\original_linearized.mat'],'X_store','tau_SRH_store');
%% Make the E-k curves 
%Given our best fits, now we make some assumptions and create the E-k
%curves at each temperature. 
%We need to load our data first
load([processing_directory '\meas_info.mat']); 
load([processing_directory '\best_fits.mat']);
label = zeros(1,1); 
defect1 = figure;
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
defect2 = figure;
tau_defect1 = figure;
tau_defect2 = figure; 
% defect3 = figure;
%I want these to be in order of increasing temperature
% for i = 1:length(info)
%     T(i) = info(i).temperature;
% end
% %Sort 
% [T,for_iteration] = sort(T); 
for i = 1:length(info)
%     index = for_iteration(i);
    index = i; 
%     best_fit = best_fits(index).one_defect; 
    best_fit = best_fits(index).two_defects;
%     best_fit = best_fits(index).three_defects;
%     Let's sort this to try to match up defects between temperatures
%     [slopes,IX] = sort(best_fit(:,1));
%     best_fit(:,1) = best_fit(IX,1); 
%     best_fit(:,2) = best_fit(IX,2); 
    [m,n] = size(best_fit);
    for j = 1:m
        [Et{index,j},k{index,j},alphanN{index,j}]=generate_Ek(best_fit(j,:),info(index).temperature+273.15,info(index).doping,type);
    end
    figure(defect1); 
    h1(i)=plot(Et{index,1},k{index,1},'-','LineWidth',2,'Color',co{i}); 
    label(i,1) = info(index).temperature; 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{index,1},1./alphanN{index,1},'-','LineWidth',2,'Color',co{i}); 
    hold all;
    figure(defect2);
    h2(i)=plot(Et{index,2},k{index,2},'-','LineWidth',2,'Color',co{i});
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{index,2},1./alphanN{index,2},'-','LineWidth',2,'Color',co{i}); 
    hold all; 
%     figure(defect3);
%     h3(i)=plot(Et{index,3},k{index,3},'-','LineWidth',2,'Color',co{i});
%     hold all;
end
%The Et vector should be the same for each defect
Et_vector = Et{index,j}; 
for i = 1:length(Et_vector)
    [temps,defects] = size(k); 
    for j = 1:defects
        for q = 1:temps
            k_now = k{q,j};
            if isnan(k_now(i))==1
               all_T(j,q) = NaN; 
            else
                all_T(j,q) = k_now(i);
            end
        end
        std_dev(i,j) = std(all_T(j,:)); 
        %We should also track the mean
        average(i,j) = mean(all_T(j,:));
    end
end
figure(defect1); 
%Now we find the minimum of standard deviation curve
% min_std = min(std_dev(:,1)); 
% Et_min = Et_vector(find(std_dev(:,1)==min_std));
% average_k = average(find(std_dev(:,1)==min_std),1); 
% hold all;
% plot([Et_min Et_min], [-100 100],'k--','LineWidth',2);
% hold all;
% plot([0 1.124], [average_k average_k],'k--','LineWidth',2);
% axis([0 1.124 0 100]);
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('Defect 1','FontSize',30); 
figure(defect2); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h2,num2str(label));
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h3,num2str(label));
title('Defect 1','FontSize',30);
figure(tau_defect2); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
xlabel('E_t-E_i [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h4,num2str(label));
title('Defect 2','FontSize',30);
% figure(defect3); 
% xlabel('E_t-E_v [eV]','FontSize',20); 
% ylabel('k [-]','FontSize',20);
% legend(h3,num2str(label));
% title('Defect 3','FontSize',30); 
figure; 
plot(Et_vector,std_dev(:,1),'LineWidth',3);
xlabel('E_t-E_v [eV]','FontSize',20);
ylabel('Standard deviation in k','FontSize',20); 
title('Defect 1','FontSize',30); 
figure; 
plot(Et_vector,std_dev(:,2),'LineWidth',3); 
xlabel('E_t-E_v [eV]','FontSize',20);
ylabel('Standard deviation in k','FontSize',20); 
title('Defect 2','FontSize',30); 
%% Plot all the E-k curves on the same plot because we don't know which defect belongs where
%We need to load our data first
load([processing_directory '\meas_info.mat']); 
load([processing_directory '\best_fits.mat']);
all_defects = figure;
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
tau_defects = figure;
defectcount = 1;
label = zeros(length(best_fits),1);
for i = 1:length(best_fits)
    index = i;
%     best_fit = best_fits(index).one_defect;
    best_fit = best_fits(index).two_defects;
%     best_fit = best_fits(index).three_defects; 
    [m,n] = size(best_fit);
    for j = 1:m
        [Et{index,j},k{index,j},alphanN{index,j}]=generate_Ek(best_fit(j,:),info(index).temperature+273.15,info(index).doping,type);
    end
    figure(all_defects); 
    h1(defectcount)=plot(Et{index,1},k{index,1},'-','LineWidth',2,'Color',co{i}); 
    label(defectcount,1) = info(index).temperature; 
    hold all; 
    figure(tau_defects); 
    h2(defectcount)=plot(Et{index,1},1./alphanN{index,1},'-','LineWidth',2,'Color',co{i}); 
    hold all;
    defectcount = defectcount+1;
    figure(all_defects);
    plot(Et{index,2},k{index,2},'--','LineWidth',2,'Color',co{i});
    hold all;
    figure(tau_defects); 
    plot(Et{index,2},1./alphanN{index,2},'--','LineWidth',2,'Color',co{i}); 
    hold all; 
%     defectcount = defectcount+1;
%     figure(all_defects);
%     plot(Et{index,3},k{index,3},'-.','LineWidth',2,'Color',co{i});
%     hold all;
%     figure(tau_defects); 
%     plot(Et{index,3},1./alphanN{index,3},'-.','LineWidth',2,'Color',co{i}); 
%     defectcount = defectcount+1;
end
figure(all_defects); 
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('All defects','FontSize',30); 
figure(tau_defects); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h2,num2str(label));
title('All defects','FontSize',30); 
%% Make best fits from Excel 
filename = [processing_directory '\IDLS two and three curve fitting.xlsm'];
load([processing_directory '\meas_info.mat']); 
fits = xlsread(filename,'Summary','B3:E11'); %Note - they should be sorted in increasing temperature
[m,n] = size(fits); 
load([processing_directory '\best_fits.mat']); %We will load this and then replace two defects fits
% for i = 1:length(info)
%     T(i) = info(i).temperature;
% end
% %Sort 
% [T,for_iteration] = sort(T); 
for i = 1:m
%     index = for_iteration(i);
    index = i;
    defect_1 = fits(i,1:2); 
    defect_2 = fits(i,3:4); 
    best_fits(index).two_defects = [defect_1;defect_2];
%     two_defects{i,1} = [defect_1;defect_2];
end
% best_fits = struct('two_defects',two_defects);
save([processing_directory '\best_fits.mat'],'best_fits');

%% %% Load the data and ask for where to crop the SRH data based on the contributions
%Load the data for processing
load([processing_directory '\lifetime_breakdown.mat']);
load([processing_directory '\meas_info.mat']);
for i = 1:length(lifetime_breakdown);
    deltan_rev = lifetime_breakdown(i).deltan; 
    tau_SRH = lifetime_breakdown(i).tau_SRH; 
    figure;
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    loglog(deltan_rev,tau_SRH.*1e6,'LineWidth',2);
    xlabel('Excess carrier density (cm^-^3)','FontSize',20);
    ylabel('Lifetime (\mus)','FontSize',20);   
    title(['Temperature = ' num2str(info(i).temperature)]);
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_highinj(deltan_rev,tau_SRH,cutoff);
    %We might always want to remove some low injection data
    disp('Select the region for cutting off the LOW injection data');
    [cutoff,nothing]=ginput(1);
    [deltan_rev,tau_SRH_rev] = remove_lowinj(deltan_rev,tau_SRH_rev,cutoff);
     %Let's store everything now
    deltan_store{i,1} = deltan_rev;
    tau_SRH_store{i,1} = tau_SRH_rev;
end
lifetime_breakdown = struct('deltan',deltan_store,'tau_SRH',tau_SRH_store);
save([processing_directory '\lifetime_breakdown.mat'],'lifetime_breakdown');
 
    
