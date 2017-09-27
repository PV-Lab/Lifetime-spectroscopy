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

%Plot the degradation curves of all samples
clc;clear all; close all; 
measurement_log = 'C:\Users\Mallory Jensen\Documents\LeTID\XRF\SERIS deg\measurement_summary.xlsx';
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\XRF\SERIS deg\Summary';
%Define the samples as they are listed in the filenames
samples = {'SALb_mid_broken','SAHb_mid','PSHb_mid','PSLb_mid'};
%Read in the data with the times
filename_details = cell(length(samples),2); 
for i = 1:length(samples)
    [num,txt,raw] = xlsread(measurement_log,samples{i}); 
    filename_details{i,1} = txt(2:end,1:2); 
    filename_details{i,2} = num; 
end
filename_start = [directory '\']; 
filename_end = {'.xlsm' '.xlsm' '.xlsm' '.xlsm'}; 

%Time intervals
% times = [0, 10:10:100, 200:100:1000, 2000:1000:10000 20000:10000:50000];

colors = {'r','g','b','m','c','y'};

zero_filenames = {[directory '\SALb_mid_broken_0s.xlsm'],...
    [directory '\SAHb_mid_0s.xlsm'],...
    [directory '\PSHb_mid_0s.xlsm'],...
    [directory '\PSLb_mid_0s.xlsm']};
%First load all of the data. Capture the injection-dependent lifetime
%curves
for i = 1:length(samples)
    figure; 
    count = 1; 
    h = []; 
    label = []; 
    %Get the number of files
    times_str = filename_details{i,1}; 
    times_str = times_str(:,2); 
    times_num = filename_details{i,2}; 
    times_num = times_num(:,4); 
    for j = 1:length(times_str)
        if times_num(j)==0
            filename = zero_filenames{i};
        else
            if i == 1 || i == 2
                 %Make the filename the way we expect it
                filename = [filename_start samples{i} '_' times_str{j} filename_end{i}];
            else
                %Make the filename the way we expect it
                filename = [filename_start samples{i} '_' times_str{j} filename_end{i}];
            end
        end
        %For the new spreadsheet
        data = xlsread(filename,'RawData','E4:G124');
        tau = data(:,1); %seconds
        deltan = data(:,3); %cm^-3
        %We need to check and make sure no deltan (x-values) are repeated. This
        %will interfere with the interpolation. 
        deltan = flipud(deltan);
        tau = flipud(tau); 
        [deltan,IX] = sort(deltan); 
        tau = tau(IX); 
        store_index = []; 
        count_index = 1; 
        for k = 1:length(deltan)
            index = find(deltan==deltan(k)); 
            [m,n] = size(index); 
            if m>1
                store_index(count_index:count_index+m-2) = index(2:end); 
                count_index = count_index+1; 
            end
        end
        deltan(store_index) = [];
        tau(store_index) = [];
        dataSave{i,j} = [deltan,tau];
        %Plot the data as we're going. We'll plot everything but we'll only
        %label the curves on the magnitude change
        if round(log10(times_num(j)))==log10(times_num(j))
            h(count) = loglog(deltan,tau.*1e6,[colors{count} '--'],'LineWidth',3);
            label(count) = times_num(j); 
            count = count+1;  
        else
            loglog(deltan,tau.*1e6,'.','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',6);
        end
        hold on; 
    end
    %Label the figure
    xlabel('Excess carrier density [cm^-^3]','FontSize',30); 
    ylabel('Lifetime [\mus]','FontSize',30); 
    title(['Sample ' samples{i}],'FontSize',30);
    legend(h,num2str(label')); 
end

%% Process and save the data after it's been loaded

%Choose one injection level and plot the lifetime at that injection level
%Now that we have data stored, it's easy enough just to cycle through it
% injection = 8e14; %cm-3
injection = 5e14; 
figure; 
for i = 1:length(samples)
    %Get the number of files
    times_str = filename_details{i,1}; 
    times_str = times_str(:,2); 
    times_num = filename_details{i,2}; 
    times_num = times_num(:,4); 
    times_num(1) = 10; 
    for j = 1:length(times_num)
        data_now = dataSave{i,j}; 
        lifetime_deg(i,j) = interp1(data_now(:,1),data_now(:,2),injection); 
    end
    loglog(times_num,lifetime_deg(i,1:length(times_num)).*1e6,'-','MarkerSize',12,'LineWidth',4); 
    hold all; 
end
xlabel('degradation time [s]','FontSize',30); 
ylabel('lifetime [\mus]','FontSize',30);
legend('low firing','high firing','PDG high firing','PDG low firing'); 
axis([1e1 1e6 3 1e2])

%Normalized degradation
figure; 
for i = 1:length(samples); 
    %Get the number of files
    times_str = filename_details{i,1}; 
    times_str = times_str(:,2); 
    times_num = filename_details{i,2}; 
    times_num = times_num(:,4); 
%     lifetime_deg_norm(i,:) = lifetime_deg(i,:)./max(lifetime_deg(i,:)); 
    lifetime_deg_norm(i,:) = lifetime_deg(i,:)./lifetime_deg(i,1);
    semilogx(times_num,lifetime_deg_norm(i,1:length(times_num)),'o','MarkerSize',12,'LineWidth',2); 
    hold all;
end
xlabel('Degradation time [s]','FontSize',30); 
ylabel('Normalized lifetime [-]','FontSize',30);
legend(samples); 

%Save the data 
save([directory '\processed_data_5e14_20170927.mat'],'dataSave','lifetime_deg_norm','lifetime_deg','filename_details','samples','times_num');
%% Plot literature degradation curves on top of an open figure
dirname = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Literature figures';
%Load Bredemeier published data
[Bred_650_1_t,Bred_650_1_tau] = read_digitized([dirname '\Bredemeier decay 1_650C.dat']);
[Bred_650_2_t,Bred_650_2_tau] = read_digitized([dirname '\Bredemeier decay 2_650C_attempt2.dat']);
[Bred_900_1_t,Bred_900_1_tau] = read_digitized([dirname '\Bredemeier decay 1_900C.dat']);
[Bred_900_2_t,Bred_900_2_tau] = read_digitized([dirname '\Bredemeier decay 2_900C_attempt2.dat']);
%Only Bredemeier 2 is relevant (AlOx/SiN). Times are in hours. Lifetimes
%are in microseconds. 
Bred_650_2_t = Bred_650_2_t.*60.*60; %s
Bred_900_2_t = Bred_900_2_t.*60.*60; %s
%Modify starting value to 1s
Bred_650_2_t(1) = 1; 
Bred_900_2_t(1) = 1; 
%Raw data
% hold all; 
% plot(Bred_650_2_t,Bred_650_2_tau); 
% hold all; 
% plot(Bred_900_2_t,Bred_900_2_tau); 
%Normalized data
Bred_650_2_tau = Bred_650_2_tau./Bred_650_2_tau(1); 
Bred_900_2_tau = Bred_900_2_tau./Bred_900_2_tau(1); 
hold all; 
plot(Bred_650_2_t,Bred_650_2_tau); 
hold all; 
plot(Bred_900_2_t,Bred_900_2_tau); 



%% Given the loaded data, try now to analyze the evolution along the degradation curve
clc;clear all; close all; 
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\rev resistivity, OC'; 
% directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\Degradation 2';
% load([directory '\processed_data.mat']);
% load([directory '\processed_data_1-1_20161103.mat']);
load([directory '\processed_data_8e14_20161122.mat']);
% sample_index = 4; 
sample_index = 2; 
%Find the fully degraded state which will correspond to the minimum
%normalized degraded lifetime
max_deg_index = find(lifetime_deg_norm(sample_index,:)==min(lifetime_deg_norm(sample_index,:))); 
%For now, max is the last index for this sample
% max_deg_index = 28; 
data_maxdeg = dataSave{sample_index,max_deg_index}; 
%Find the fully UNdegraded state which is the state that we started in
data_mindeg = dataSave{sample_index,1}; 
%Plot the two lifetimes together
figure;
loglog(data_mindeg(:,1),data_mindeg(:,2),'LineWidth',3); 
hold all;
loglog(data_maxdeg(:,1),data_maxdeg(:,2),'LineWidth',3);
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
legend('Fully undegraded','Fully degraded'); 

%Try calculating the SRH lifetime always relative to the initial state

% to_calc = [80 85 87 91 92 93 96]; 
% to_calc = [20 29 39 43 46 47];
% labels = [500 1000 10000 100000 200021 300297 380757]; %degradation
% to_calc = [6 11 20 24]; %deg 2
% to_calc = [8 11 16 20 24 28]; %deg semifabricate
[num_samples,num_measurements] = size(dataSave); 
times_num = filename_details{sample_index,2}; times_num = times_num(:,2); 
SRHfig = figure;
taufig = figure;
h2(1)=loglog(data_mindeg(:,1),data_mindeg(:,2),'LineWidth',3); 
hold all;
for i = 1:length(to_calc)
    datanow = dataSave{sample_index,to_calc(i)}; 
    %Plot the raw lifetime for publication
    figure(taufig); 
    h2(i+1) = loglog(datanow(:,1),datanow(:,2),'LineWidth',3); 
    hold all; 
    %We need to interpolate the lifetime
    tau_measure = interp1(datanow(:,1),datanow(:,2),data_mindeg(:,1)); 
    tau_SRH{i} = ((1./tau_measure)-(1./data_mindeg(:,2))).^(-1);
    figure(SRHfig)
    h(i)=loglog(data_mindeg(:,1),tau_SRH{i},'LineWidth',3); 
    labels(i) = times_num(to_calc(i)); 
    hold all;
end
figure(SRHfig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('SRH Lifetime [s]','FontSize',20);
legend(h,num2str(labels')); 
figure(taufig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
legend(h,num2str([0;labels'])); 

%Now let's take this SRH lifetime and try fitting it!
doping = 6.9e15; 
T = 300; 
type = 'p'; 
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,doping,type); 
fit_tries = 1e6; 
for i = 1:length(to_calc)
    tau= tau_SRH{i}; 
    deltan = data_mindeg(:,1); 
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
     if type == 'p'
        X = (n0+deltan_rev)./(p0+deltan_rev);
    elseif type == 'n'
        X = (p0+deltan_rev)./(n0+deltan_rev);
     end
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_rev.*1e6,25,directory,fit_tries);
end
%Select two defects and generate the E_k curves
defect1 = figure;
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
defect2 = figure;
tau_defect1 = figure;
tau_defect2 = figure; 
for i = 1:length(to_calc)
    best_fit = two_defects{i,1};
    [slopes,IX] = sort(best_fit(:,1));
    best_fit(:,1) = best_fit(IX,1); 
    best_fit(:,2) = best_fit(IX,2); 
    for j = 1:length(best_fit)
        [Et{i,j},k{i,j},alphanN{i,j}]=generate_Ek(best_fit(j,:),T,doping,type);
    end
    figure(defect1); 
    h1(i)=plot(Et{i,1},k{i,1},'-','LineWidth',2,'Color',co{i}); 
    label(i,1) = T; 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',2,'Color',co{i}); 
    hold all;
    figure(defect2);
    h2(i)=plot(Et{i,2},k{i,2},'-','LineWidth',2,'Color',co{i});
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{i,2},1./alphanN{i,2},'-','LineWidth',2,'Color',co{i}); 
    hold all; 
end
figure(defect1); 
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('Defect 1','FontSize',30); 
figure(defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h2,num2str(label));
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [\mus]','FontSize',20);
legend(h3,num2str(label));
title('Defect 1','FontSize',30);
figure(tau_defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [\mus]','FontSize',20);
legend(h4,num2str(label));
title('Defect 2','FontSize',30);
save([directory '\fitted_defect_parameters_addedRegen.mat'],'Et','alphanN','k','two_defects','T','doping','type','samples','labels');
%% Vary the lifetime assess the uncertainty in k values
clc;clear all; close all; 
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC';
load([directory '\processed_data.mat']);
%Find the fully degraded state which will correspond to the minimum
%normalized degraded lifetime
max_deg_index = find(lifetime_deg_norm(2,:)==min(lifetime_deg_norm(2,:))); 
data_maxdeg = dataSave{2,max_deg_index}; 
%Find the fully UNdegraded state which is the state that we started in
data_mindeg = dataSave{2,1}; 
%Plot the two lifetimes together
figure;
loglog(data_mindeg(:,1),data_mindeg(:,2),'LineWidth',3); 
hold all;
loglog(data_maxdeg(:,1),data_maxdeg(:,2),'LineWidth',3);
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
legend('Fully undegraded','Fully degraded'); 

%Try calculating the SRH lifetime always relative to the initial state
to_calc = [20 29 39 43 46 47];
error = [-18 -9 -7 -7 -7 -7];
[num_samples,num_measurements] = size(dataSave); 
times_num = filename_details{2,2}; times_num = times_num(:,4); 
SRHfig = figure;
taufig = figure;
h2(1)=loglog(data_mindeg(:,1),data_mindeg(:,2),'LineWidth',3); 
hold all;
for i = 1:length(to_calc)
    datanow = dataSave{2,to_calc(i)}; 
    %Plot the raw lifetime for publication
    figure(taufig); 
    h2(i+1) = loglog(datanow(:,1),datanow(:,2),'LineWidth',3); 
    hold all; 
    %We need to interpolate the lifetime
    tau_measure = interp1(datanow(:,1),datanow(:,2),data_mindeg(:,1)); 
    tau_SRH_raw = ((1./tau_measure)-(1./data_mindeg(:,2))).^(-1);
    tau_SRH{i} = tau_SRH_raw.*(1+(error(i)/100)); 
    figure(SRHfig)
    h(i)=loglog(data_mindeg(:,1),tau_SRH{i},'LineWidth',3); 
    labels(i) = times_num(to_calc(i)); 
    hold all;
end
figure(SRHfig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('SRH Lifetime [s]','FontSize',20);
legend(h,num2str(labels')); 
figure(taufig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
legend(h,num2str([0;labels'])); 

%Now let's take this SRH lifetime and try fitting it!
doping = 6.9e15; 
T = 300; 
type = 'p'; 
[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,doping,type); 
fit_tries = 1e6; 
for i = 1:length(to_calc)
    tau= tau_SRH{i}; 
    deltan = data_mindeg(:,1); 
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
     if type == 'p'
        X = (n0+deltan_rev)./(p0+deltan_rev);
    elseif type == 'n'
        X = (p0+deltan_rev)./(n0+deltan_rev);
     end
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_rev.*1e6,25,directory,fit_tries);
end
%Select two defects and generate the E_k curves
defect1 = figure;
co={[0 0 0]; [0.5 0 0.9]; [0 0 1]; [0 1 1]; [0 1 0];  [1 1 0]; [1 0.6 0]; [1 0 0]; [0.8 0.5 0]};
defect2 = figure;
tau_defect1 = figure;
tau_defect2 = figure; 
for i = 1:length(to_calc)
    best_fit = two_defects{i,1};
    [slopes,IX] = sort(best_fit(:,1));
    best_fit(:,1) = best_fit(IX,1); 
    best_fit(:,2) = best_fit(IX,2); 
    for j = 1:length(best_fit)
        [Et{i,j},k{i,j},alphanN{i,j}]=generate_Ek(best_fit(j,:),T,doping,type);
    end
    figure(defect1); 
    h1(i)=plot(Et{i,1},k{i,1},'-','LineWidth',2,'Color',co{i}); 
    label(i,1) = T; 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',2,'Color',co{i}); 
    hold all;
    figure(defect2);
    h2(i)=plot(Et{i,2},k{i,2},'-','LineWidth',2,'Color',co{i});
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{i,2},1./alphanN{i,2},'-','LineWidth',2,'Color',co{i}); 
    hold all; 
end
figure(defect1); 
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h1,num2str(label));
title('Defect 1','FontSize',30); 
figure(defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('k [-]','FontSize',20);
legend(h2,num2str(label));
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [\mus]','FontSize',20);
legend(h3,num2str(label));
title('Defect 1','FontSize',30);
figure(tau_defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [\mus]','FontSize',20);
legend(h4,num2str(label));
title('Defect 2','FontSize',30);
save([directory '\fitted_defect_parameters_werror_neg.mat'],'Et','alphanN','k','two_defects','T','doping','type','samples','labels');



%% Try some other manipulations of the data for presentation purposes
clear all; close all; 
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC';
load([directory '\fitted_defect_parameters.mat']); 
%Now plot the k-values and tau_n0 values for each defect at midgap
[num_meas,num_defects] = size(k); 
for i = 1:num_defects
    midgap = 1.124/2; 
    tau_n0_fig = figure;
    k_fig = figure;
    for j = 1:num_meas
        Et_now = Et{j,i};
        mg_index = find(min(abs(Et_now-midgap))==abs(Et_now-midgap)); %closest index to midgap
        k_now = k{j,i}; 
        alphanN_now = alphanN{j,i};
        figure(k_fig); 
        plot(labels(j),k_now(mg_index),'x','MarkerSize',10,'LineWidth',2); 
        hold on;
        figure(tau_n0_fig); 
        plot(labels(j),alphanN_now(mg_index),'x','MarkerSize',10,'LineWidth',2); 
        hold on;
    end
    figure(k_fig); 
    xlabel('Degradation time [s]','FontSize',30); 
    ylabel('\sigma_n/\sigma_p [-]','FontSize',30); 
    title(['Defect ' num2str(i)],'FontSize',30); 
    figure(tau_n0_fig); 
    xlabel('Degradation time [s]','FontSize',30); 
    ylabel('N_t\sigma_nv_th','FontSize',30); 
    title(['Defect ' num2str(i)],'FontSize',30); 
end

%% Post-processing - midgap defect values with error bars
clear all; close all; clc;
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\Processed degraded';
load([directory '\fitted_defect_parameters_deg.mat']); 
labels = [500 1000 10000 100000 200021 300297 380757]; %degradation
%Now plot the k-values and tau_n0 values for each defect at midgap
[num_meas,num_defects] = size(k); 
alphanN_store = zeros(num_meas,3); %for Tidd, Mo, W
% for i = 1:num_defects
for i = 2 %we only care about defect 2
    midgap = 1.124/2; 
    tau_n0_fig = figure;
    k_fig = figure;
    for j = 1:num_meas
        Et_now = Et{j,i};
        mg_index = find(min(abs(Et_now-midgap))==abs(Et_now-midgap)); %closest index to midgap
        k_now = k{j,i}; 
        alphanN_now = alphanN{j,i};
        k_store(j,i) = k_now(mg_index);
%         alphanN_store(j,i) = alphanN_now(mg_index);
        alphanN_store(j,1) = alphanN_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        alphanN_store(j,2) = alphanN_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        alphanN_store(j,3) = alphanN_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
    end
end
%I only care about defect 2 for the moment
k_def2 = k_store(:,2); 
% alphanN_def2 = alphanN_store(:,2); 
alphanN_def2 = alphanN_store; 
%For k, we define a general error based on previous experiments
k_error = k_def2.*0.16; %+/- 16%
figure(k_fig); 
colors={[0.87 0.49 0]; [0.75 0.75 0]; [0 0.5 0]; [0 0.75 0.75]; [0 0 1];  [0.75 0 0.75]; [0 0 0]};%degradation
% colors={[1 0 0]; [0.87 0.49 0]; [0.75 0.75 0]; [0 0.5 0]; [0 0.75 0.75]; [0 0 1];  [0.75 0 0.75]; [0.6 0.59 0]};%regeneration

% colors = {[0 0.4980 0],[1 0 0],[0 0.749 0.749],[0.749 0 0.749],[0.749 0.749 0],[0.2471 0.2471 0.2471]};
for i = 1:length(k_def2)
    errorbar(labels(i)',k_def2(i),k_error(i),'x','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
%     errorbar(labels(i)',k_def2(i),k_error(i),'x','MarkerSize',12,'LineWidth',3); 
%     errorbar(labels(i),k_def2(i),k_error(i),'x','MarkerSize',12,'LineWidth',3); 
    hold on;
end
%Plot bounds from previous
hold on; 
semilogx([1 1e6],[26 26],'--','LineWidth',2,'Color',[0.5 0.5 0.5]); 
hold on; 
semilogx([1 1e6], [36 36],'--','LineWidth',2,'Color',[0.5 0.5 0.5]); 
xlabel('Degradation time [s]','FontSize',30); 
ylabel('k [-]','FontSize',30); 
axis([0 1e6 0 50]);
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
%For tau_n0, we can define the error based on perturbations in the measured
%lifetime
% error_upper = [1.13 1.09 1.07 1.06 1.06 1.06];
% error_lower = [0.82 0.9 0.93 0.93 0.93 0.93]; 
error_upper = 1.15.*ones(1,num_meas); 
error_lower = 0.85.*ones(1,num_meas); 
%Let's try plotting for different possible defects
sigman_1 = 1.5e-15; %Tidd
sigman_2 = 1.6e-14; %Mo
sigman_3 = 1.7e-14; %W
%Modify the units so that they are in seconds rather than microseconds
alphanN_def2 = alphanN_def2.*1e6; %1/s
%Let's take out the thermal velocity for clarity. These terms are
%Nt*sigma_n*vth_e
vth_e = 2.05e7; %cm/s
alphanN_def2 = alphanN_def2./vth_e; %Nt*sigma_n, units cm-1
%Now also take out the different sigma values
anN_1 = alphanN_def2(:,1)./sigman_1; 
anN_2 = alphanN_def2(:,2)./sigman_2; 
anN_3 = alphanN_def2(:,3)./sigman_3; 
%Now we find the error founds
U_1 = anN_1.*error_upper'; %1/cm3
L_1 = anN_1.*error_lower'; %1/cm3
U_2 = anN_2.*error_upper'; %1/cm3
L_2 = anN_2.*error_lower'; %1/cm3
U_3 = anN_3.*error_upper'; %1/cm3
L_3 = anN_3.*error_lower'; %1/cm3
% %Modify the units so that they are in seconds rather than microseconds
% U = U.*1e6; %1/s
% L = L.*1e6; %1/s
% alphanN_def2.*1e6; %1/s
% %Let's take out the thermal velocity for clarity. These terms are
% %Nt*sigma_n*vth_e
% vth_e = 2.05e7; %cm/s
% U = U./vth_e; %Nt*sigma_n, units cm-1
% L = L./vth_e; 
% alphanN_def2 = alphanN_def2./vth_e; 
figure(tau_n0_fig); 
for i = 1:length(anN_1)
    h(1)=loglog(labels(i),anN_1(i),'s','MarkerSize',12,'LineWidth',3,'Color',colors{i});
    hold all; 
    h(2) = loglog(labels(i),anN_2(i),'o','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
    hold all;
    h(3)=loglog(labels(i),anN_3(i),'x','MarkerSize',14,'LineWidth',3,'Color',colors{i}); 
    hold all; 
%     h(1)=loglog(labels(i),anN_1(i),'x','MarkerSize',12,'LineWidth',2,'Color',colors{i}); 
%     hold all;
%     h(2)=loglog(labels(i),anN_3(i),'s','MarkerSize',12,'LineWidth',2,'Color',colors{i}); 
%     h(1)=errorbar(labels(i),anN_1(i),L_1(i),U_1(i),'x','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
%     hold on;
%     errorbar(labels(i),anN_2(i),L_2(i),U_2(i),'o','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
%     hold all;
%     h(2)=errorbar(labels(i),anN_3(i),L_3(i),U_3(i),'s','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
end
legend(h,'Ti','Mo','W');
set(get(h(1),'Parent'),'YScale','log');
xlabel('Degradation time [s]','FontSize',30); 
ylabel('N_t [cm^{-3}]','FontSize',30); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
        
%% Process PLI
clear all; close all; clc;
%Define delimiter and header in order to read txt files
delimiterIn = ',';
headerlinesIn = 0;
exposure = [10 30]; 
samples = {'64-5' '69-5'};
LP = 50; 
filename_start = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\PLI';
filename_end = 'LP_1.txt';
times = [981787 1046047 1105807 1136527 1177510 1207810 1258810 1340770 1377370 1422190 1479390 1530090 1563810]; 
folders = {'981787s' '1046047s' '1105807s' '1136527s' '1177510s' '1207810s' '1258810s' '1340770s' '1377370s' '1422190s' '1479390s' '1530090s' '1563810s'}; 
ref_filename = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\PLI\reference\68-5_10s_50LP_1.txt';
PLref = importdata(ref_filename,delimiterIn,headerlinesIn);
PLref = PLref(:,2:end)./10; %counts/second
[PLref,x_crop,y_crop]=crop_PL(PLref);
minvalue_ref = min(min(PLref));
maxvalue_ref = max(max(PLref));
scaled_ref = (PLref-minvalue_ref)./(maxvalue_ref-minvalue_ref);
h=figure; 
%we need to plot the reference file
imagesc(scaled_ref);
axis('image');
colormap('gray');
caxis([0 1]);
colorbar;
title('Initial'); 
axis off; 
print(h,'-dpng','-r0',[filename_start '\Initial_PLI.png']);
for i = 1:length(samples)
    for j = 1:length(times); 
        %Make filename
        filename = [filename_start '\' folders{j} '\' samples{i} '_' num2str(times(j)) 's_' num2str(exposure(i)) 's_' num2str(LP) filename_end]; 
        %Load data
        PLmap = importdata(filename,delimiterIn,headerlinesIn);
        PLmap = PLmap(:,2:end)./exposure(i); %counts/second
        [PLmap,x_crop,y_crop]=crop_PL(PLmap);
        minvalue = min(min(PLmap));
        maxvalue = max(max(PLmap));
        scaled = (PLmap-minvalue)./(maxvalue-minvalue);
        h=figure; 
        imagesc(scaled);
        axis('image');
        colormap('gray');
        caxis([0 1]);
        colorbar;
        title([samples{i} ', ' folders{j}]); 
        axis off; 
        print(h,'-dpng','-r0',[filename_start '\' samples{i} '_' folders{j} '.png']);
    end
end
%% Assess lifetime fitting with same parameters as previous publication
clc;clear all; close all; 
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\rev resistivity, OC';
% load([directory '\processed_data_1-1_20161103.mat']);
load([directory '\processed_data_8e14_20161122.mat']);
sample_index = 2;
%Find the fully degraded state which will correspond to the minimum
%normalized degraded lifetime
max_deg_index = find(lifetime_deg_norm(sample_index,:)==min(lifetime_deg_norm(sample_index,:))); 
% max_deg_index = 22; 
data_maxdeg = dataSave{sample_index,max_deg_index}; 
%Find the fully UNdegraded state which is the state that we started in
data_mindeg = dataSave{sample_index,1}; 
%Plot the two lifetimes together
figure;
loglog(data_mindeg(:,1),data_mindeg(:,2),'LineWidth',3); 
hold all;
loglog(data_maxdeg(:,1),data_maxdeg(:,2),'LineWidth',3);
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
legend('Fully undegraded','Fully degraded'); 

%Try calculating the SRH lifetime always relative to the initial state
%These correspond to these times:

%DEGRADATION
%15=500s or 8:20, 20=1000s or 16:40, 
%29=10000s or 2:46:40, 39=100000s or 24:224:160, 43=200021s or 51:270:221,
%46=300297 or 79:261:237, 47=380757 or 101:282:237
% to_calc = [15 20 29 39 43 46 47];
to_calc = [12 13 14];

%REGENERATION
% to_calc = [47 55 61 70 74 75 77 79]; 
% to_calc = [80 85 87 91 92 93 96]; %added values for review
% to_calc = [80 87 91 92 93 96]; %added values for review
% to_calc = [80 91 92 93 96]; %added values for review

%DEGRADATION 2
% to_calc = [22];

%SEMIFABRICATE DEGRADATION
% to_calc = [12 13 14 15 16 20 24 28]; %deg semifabricate

[num_samples,num_measurements] = size(dataSave); 
times_num = filename_details{2,2}; times_num = times_num(:,4); 
% times_num = filename_details{sample_index,2}; times_num = times_num(:,1); 
SRHfig = figure;
taufig = figure;
h2(1)=loglog(data_mindeg(:,1),data_mindeg(:,2),'LineWidth',3); 
hold all;
for i = 1:length(to_calc)
    datanow = dataSave{sample_index,to_calc(i)}; 
    %Plot the raw lifetime for publication
    figure(taufig); 
    h2(i+1) = loglog(datanow(:,1),datanow(:,2),'LineWidth',3); 
    hold all; 
    %We need to interpolate the lifetime
    tau_measure = interp1(datanow(:,1),datanow(:,2),data_mindeg(:,1)); 
    tau_SRH{i} = ((1./tau_measure)-(1./data_mindeg(:,2))).^(-1);
    figure(SRHfig)
    h(i)=loglog(data_mindeg(:,1),tau_SRH{i},'LineWidth',3); 
    labels(i) = times_num(to_calc(i)); 
    hold all;
end
figure(SRHfig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('SRH Lifetime [s]','FontSize',20);
title('Regeneration','FontSize',20); 
legend(h,num2str(labels')); 
% print(SRHfig,'-dpng','-r0',[directory '\SRH_lifetimes_deg_semifab.png']);
% hgsave(SRHfig,[directory '\SRH_lifetimes_deg_semifab']);
figure(taufig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
title('Regeneration','Fontsize',20);
legend(h2,num2str([0;labels'])); 
% print(taufig,'-dpng','-r0',[directory '\Measured_lifetimes_deg_semifab.png']);
% hgsave(taufig,[directory '\Measured_lifetimes_deg_semifab']);
%% Continue with the fitting process after running previous section

%Now let's take this SRH lifetime and try fitting it!
doping = 9.09e15; 
T = 300; 
type = 'p'; 
%hard code the parameters as used for publication
n0 = 1.02e4; 
p0 =  9.09e15; 
NC = 3e19; 
NV = 1e19; 
vth_e = 2.05e7; 
vth_h = 1.69e7; 
Eg = 1.1242; 
k_B = 8.62e-5; 
fit_tries = 1e6; 
for i = 1:length(to_calc)
    tau= tau_SRH{i}; 
    deltan = data_mindeg(:,1); 
    %Plot the data and ask the user where the cut off in high injection
    figure;
    loglog(deltan,tau,'.');
    disp('Select the region for cutting off the HIGH injection data');
    [cutoff,nothing]=ginput(1);
%     cutoff = 1e16; 
    [deltan_rev,tau_rev] = remove_highinj(deltan,tau,cutoff);
    %We might always want to remove some low injection data
%     disp('Select the region for cutting off the LOW injection data');
%     [cutoff,nothing]=ginput(1);
    cutoff = 2e14; 
    [deltan_rev,tau_rev] = remove_lowinj(deltan_rev,tau_rev,cutoff);
    hold all;
    loglog(deltan_rev,tau_rev,'+');
    legend('Before cutoff','After cutoff'); 
     if type == 'p'
        X = (n0+deltan_rev)./(p0+deltan_rev);
    elseif type == 'n'
        X = (p0+deltan_rev)./(n0+deltan_rev);
     end
    xlswrite([directory '\Linearized_data_addedEarlyDeg.xlsx'],[X,tau_rev],['Sheet' num2str(i)]); 
    [one_defect{i,1},MSE_one{i,1},two_defects{i,1},MSE_two{i,1},three_defects{i,1},MSE_three{i,1},all_parameters_store{i,1},all_MSE_store{i,1}] = fit_murphy_master(X,tau_rev.*1e6,25,directory,fit_tries);
    to_write = zeros(6,3); 
    to_write(1:2,1) = one_defect{i,1}';
    twodef = two_defects{i,1}; 
    to_write(1:2,2) = twodef(1,:)';
    to_write(3:4,2) = twodef(2,:)';
    threedef = three_defects{i,1}; 
    to_write(1:2,3) = threedef(1,:)';
    to_write(3:4,3) =threedef(2,:)';
    to_write(5:6,3) = threedef(3,:)';
    xlswrite([directory '\Linearized_data_addedEarlyDeg.xlsx'],to_write,['Sheet' num2str(i)],'C1:E6'); 
end
%Now pause and refine fits in Excel. 
%% After fits refined in Excel, make Ek curves
clear all; close all; 
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\rev resistivity, OC\Processed degraded\added for review';
% directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\Processed degraded';
% directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\Degradation 2\Analyzing 1-1';
% to_calc = [15 20 29 39 43 46 47]; %degradation
% to_calc = [47 55 61 70 74 75 77 79]; %regeneration
% to_calc = [12 13 14 15 16 20 24 28]; %deg semifabricate
% to_calc = [12 16 20 24 28]; %deg semifabricate
% to_calc = [11 13 16 20 24 22];%Degradation 2
% to_calc = [80 85 87 91 92 93 96]; %added values for review
to_calc = [11 12 13 14]; %added early degradation values for review
% label = {'8.3 min';'16.7 min';'2.8 hr';'27.8 hr';'55.6 hr';'83.4 hr';'105.8 hr'}; %degradation
% label = {'105.8 hr';'214.9 hr';'290.6 hr';'425.0 hr';'572.5 hr';'646.0 hr';'715.2 hr';'759.4 hr'}; %regeneration
% label = {'8 min';'17 min';'3 hr';'28 hr';'56 hr';'83 hr';'106 hr'}; %degradation
% label = {'106 hr';'215 hr';'290 hr';'425 hr';'572 hr';'646 hr';'715 hr';'759 hr'}; %regeneration
% label = {'17 min';'50 min';'1.7 hr';'2.8 hr';'14.3 hr';'8.3 hr'};
% label = {'2752980s','3012600s','3210480s','3963900s','4213650s','4554330s','5002890s'};
label = {'100s','200s','300s','400s'};
% label = {'2000s';'6000s';'10000s';'51642s';'102375s';'3000s';'4000s';'5000s'};
%Now let's take this SRH lifetime and try fitting it!
doping = 9.09e15; 
T = 300; 
type = 'p'; 
%hard code the parameters as used for publication
n0 = 1.02e4; 
p0 =  9.09e15; 
NC = 3e19; 
NV = 1e19; 
vth_e = 2.05e7; 
vth_h = 1.69e7; 
Eg = 1.1242; 
k_B = 8.61733238e-5;
%Make best fits from Excel
filename = [directory '\IDLS two and three curve fitting.xlsm'];
fits = xlsread(filename,'Summary','B3:E11'); %They should be in the same order as those in the MATLAB file
[m,n] = size(fits); 
two_defects = cell(m,1); 
for i = 1:m
    index = i;
    defect_1 = fits(i,1:2); 
    defect_2 = fits(i,3:4); 
    two_defects{i,1} = [defect_1;defect_2];
end
best_fits = struct('two_defects',two_defects);
save([directory '\best_fits.mat'],'best_fits');

%Select two defects and generate the E_k curves
defect1 = figure;
% co={[222 125 0]; [191 191 0]; [0 127 0]; [0 191 191]; [0 0 255];  [191 0 191]; [153 151 0]};%degradation
% co={[255 0 0]; [222 125 0]; [191 191 0]; [0 127 0]; [0 191 191]; [0 0 255];  [191 0 191]; [153 51 0]};%regeneration
% co={[0.87 0.49 0]; [0.75 0.75 0]; [0 0.5 0]; [0 0.75 0.75]; [0 0 1];  [0.75 0 0.75]; [0 0 0]};%degradation
co={[1 0 0]; [0.87 0.49 0]; [0.75 0.75 0]; [0 0.5 0]; [0 0.75 0.75]; [0 0 1];  [0.75 0 0.75]; [0.6 0.2 0]};%regeneration
defect2 = figure;
tau_defect1 = figure;
tau_defect2 = figure; 
Et= cell(length(to_calc),2); 
k= cell(length(to_calc),2); 
alphanN= cell(length(to_calc),2); 
for i = 1:length(to_calc)
    best_fit_hold = two_defects{i,1};
%     [slopes,IX] = sort(best_fit_hold(:,1));
%     best_fit_hold(:,1) = best_fit_hold(IX,1); 
%     best_fit_hold(:,2) = best_fit_hold(IX,2); 
%     for j = 1:length(best_fit_hold)
%         [Et{i,j},k{i,j},alphanN{i,j}]=generate_Ek(best_fit_hold(j,:),T,doping,type);
%     end
    for j = 1:length(best_fit_hold)
        best_fit_now = best_fit_hold(j,:); 
        %Define the energy levels for evaluation
        Et_now = linspace(0,Eg,250); %eV
        Q = zeros(size(Et_now)); 
        alphanN_now = zeros(size(Et_now));
        k_now = zeros(size(Et_now)); 
        A = best_fit_now(1)+best_fit_now(2); %X -> 1
        B = best_fit_now(2); %X -> 0
        C = best_fit_now(1)/A; %slope/X -> 1
        for l = 1:length(Et_now)
            %Calculate n1
            n1 = NC*exp(-(Eg-Et_now(l))/(k_B*T)); 
            %Calculate p1
            p1 = NV*exp(-Et_now(l)/(k_B*T)); 
            %Calculate the Q values for these defects
            Q(l) = (C+(p1/p0))/(1-(n1/p0)-C);
            %Calculate the quantity alphan*Nt for these defects
            alphanN_now(l) = (1/B)*(1+((1/p0)*((Q(l)*n1)+p1))); 
            %Calculate the k values for these defects
            k_now(l) = Q(l)*vth_h/vth_e; 
        end
        %Get rid of any negative k values
        indices = find(k_now<0); 
        if length(indices)<length(k_now)
            Et_now(indices) = []; 
            k_now(indices) = []; 
            alphanN_now(indices) = [];  
        else
            warning(['All of the entries were negative for defect ' num2str(j) ', time ' label{i}]);
        end
        Et{i,j} = Et_now; 
        k{i,j} = k_now; 
        alphanN{i,j} = alphanN_now; 
    end
    figure(defect1); 
%     h1(i)=plot(Et{i,1},k{i,1},'-','LineWidth',3,'Color',co{i}); 
    plot(Et{i,1},k{i,1},'-','LineWidth',3,'Color',co{i}); 
%     label(i,1) = T; 
    hold all; 
    figure(tau_defect1); 
%     h3(i)=plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',3,'Color',co{i});
    plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',3,'Color',co{i});
    hold all;
    figure(defect2);
    h2(i)=plot(Et{i,2},k{i,2},'-','LineWidth',3,'Color',co{i});
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{i,2},1./alphanN{i,2},'-','LineWidth',3,'Color',co{i}); 
    hold all; 
end
figure(defect1); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('k [-]','FontSize',30);
% legend(h1,label);
title('Defect 1','FontSize',30); 
figure(defect2);
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 50]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('k [-]','FontSize',30);
legend(h2,label);
%Let's also plot the identified range from previous paper 26-36 (k)
hold all; 
plot([0 1.13],[26 26],'--','LineWidth',2,'Color',[0.5 0.5 0.5]); 
hold all; 
plot([0 1.13],[36 36],'k--','LineWidth',2,'Color',[0.5 0.5 0.5]);
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('\tau_{n0} [\mus]','FontSize',30);
% legend(h3,label);
title('Defect 1','FontSize',30);
figure(tau_defect2); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('\tau_{n0} [\mus]','FontSize',30);
legend(h4,label);
title('Defect 2','FontSize',30);
save([directory '\fitted_defect_parameters_earlyDeg.mat'],'Et','alphanN','k','two_defects','T','doping','type','label');
%% Summarize defect parameters throughout both degradation and regeneration
clear all; close all; clc;
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\Processed degraded';
load([directory '\fitted_defect_parameters_deg.mat']); 
labels = [500 1000 10000 100000 200021 300297 380757]; %degradation
%Now plot the k-values and tau_n0 values for each defect at midgap
[num_meas,num_defects] = size(k); 
alphanN_store = zeros(num_meas,3); %for Tidd, Mo, W
% for i = 1:num_defects
for i = 2 %we only care about defect 2
    midgap = 1.124/2; 
    tau_n0_fig = figure;
    k_fig = figure;
    for j = 1:num_meas
        Et_now = Et{j,i};
        mg_index = find(min(abs(Et_now-midgap))==abs(Et_now-midgap)); %closest index to midgap
        k_now = k{j,i}; 
        alphanN_now = alphanN{j,i};
%         k_store(j,i) = k_now(mg_index);
        k_store(j,1) = k_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        k_store(j,2) = k_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        k_store(j,3) = k_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
        alphanN_store(j,1) = alphanN_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        alphanN_store(j,2) = alphanN_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        alphanN_store(j,3) = alphanN_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
    end
end
%I only care about defect 2 for the moment
% k_def2 = k_store(:,2); 
k_def2 = k_store; 
% alphanN_def2 = alphanN_store(:,2); 
alphanN_def2 = alphanN_store; 
%For k, we define a general error based on previous experiments
k_error = k_def2.*0.16; %+/- 16%
figure(k_fig); 
for i = 1:length(k_def2)
    errorbar(labels(i)',k_def2(i),k_error(i),'kx','MarkerSize',12,'LineWidth',3); 
    hold on;
end
%For tau_n0, we can define the error based on perturbations in the measured
%lifetime
error_upper = 1.15.*ones(1,num_meas); 
error_lower = 0.85.*ones(1,num_meas); 
%Let's try plotting for different possible defects
sigman_1 = 1.5e-15; %Tidd
sigman_2 = 1.6e-14; %Mo
sigman_3 = 1.7e-14; %W
%Modify the units so that they are in seconds rather than microseconds
alphanN_def2 = alphanN_def2.*1e6; %1/s
%Let's take out the thermal velocity for clarity. These terms are
%Nt*sigma_n*vth_e
vth_e = 2.05e7; %cm/s
alphanN_def2 = alphanN_def2./vth_e; %Nt*sigma_n, units cm-1
%Now also take out the different sigma values
anN_1 = alphanN_def2(:,1)./sigman_1; 
anN_2 = alphanN_def2(:,2)./sigman_2; 
anN_3 = alphanN_def2(:,3)./sigman_3; 
%Now we find the error founds
U_1 = anN_1.*error_upper'; %1/cm3
L_1 = anN_1.*error_lower'; %1/cm3
U_2 = anN_2.*error_upper'; %1/cm3
L_2 = anN_2.*error_lower'; %1/cm3
U_3 = anN_3.*error_upper'; %1/cm3
L_3 = anN_3.*error_lower'; %1/cm3
figure(tau_n0_fig); 
for i = 1:length(anN_1)
    h(1)=loglog(labels(i),anN_1(i),'ks','MarkerSize',12,'LineWidth',3);
    hold all; 
    h(2) = loglog(labels(i),anN_2(i),'ko','MarkerSize',12,'LineWidth',3); 
    hold all;
    h(3)=loglog(labels(i),anN_3(i),'kx','MarkerSize',14,'LineWidth',3); 
    hold all;  
end
%Now we need to add the regeneration data in the exact same way
%to_calc = [47 55 61 70 74 75 77 79]; 
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\Processed regeneration\Excluding low injection data';
load([directory '\fitted_defect_parameters_regen.mat']); 
labels = [380757 826037 1105807 1530090 2061150 2325570 2574840 2733780]; %regeneration
%Now plot the k-values and tau_n0 values for each defect at midgap
[num_meas,num_defects] = size(k); 
alphanN_store = zeros(num_meas,3); %for Tidd, Mo, W
for i = 2 %we only care about defect 2
    midgap = 1.124/2; 
    for j = 1:num_meas
        Et_now = Et{j,i};
        mg_index = find(min(abs(Et_now-midgap))==abs(Et_now-midgap)); %closest index to midgap
        k_now = k{j,i}; 
        alphanN_now = alphanN{j,i};
        k_store(j,i) = k_now(mg_index);
        alphanN_store(j,1) = alphanN_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        alphanN_store(j,2) = alphanN_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        alphanN_store(j,3) = alphanN_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
    end
end
%I only care about defect 2 for the moment
k_def2 = k_store(:,2); 
alphanN_def2 = alphanN_store; 
%For k, we define a general error based on previous experiments
k_error = k_def2.*0.16; %+/- 16%
figure(k_fig); 
for i = 1:length(k_def2)
    errorbar(labels(i)',k_def2(i),k_error(i),'kx','MarkerSize',12,'LineWidth',3); 
    hold on;
end
%For tau_n0, we can define the error based on perturbations in the measured
%lifetime
error_upper = 1.15.*ones(1,num_meas); 
error_lower = 0.85.*ones(1,num_meas); 
%Let's try plotting for different possible defects
sigman_1 = 1.5e-15; %Tidd
sigman_2 = 1.6e-14; %Mo
sigman_3 = 1.7e-14; %W
%Modify the units so that they are in seconds rather than microseconds
alphanN_def2 = alphanN_def2.*1e6; %1/s
%Let's take out the thermal velocity for clarity. These terms are
%Nt*sigma_n*vth_e
vth_e = 2.05e7; %cm/s
alphanN_def2 = alphanN_def2./vth_e; %Nt*sigma_n, units cm-1
%Now also take out the different sigma values
anN_1 = alphanN_def2(:,1)./sigman_1; 
anN_2 = alphanN_def2(:,2)./sigman_2; 
anN_3 = alphanN_def2(:,3)./sigman_3; 
%Now we find the error founds
U_1 = anN_1.*error_upper'; %1/cm3
L_1 = anN_1.*error_lower'; %1/cm3
U_2 = anN_2.*error_upper'; %1/cm3
L_2 = anN_2.*error_lower'; %1/cm3
U_3 = anN_3.*error_upper'; %1/cm3
L_3 = anN_3.*error_lower'; %1/cm3
figure(tau_n0_fig); 
for i = 1:length(anN_1)
    h(1)=loglog(labels(i),anN_1(i),'ks','MarkerSize',12,'LineWidth',3);
    hold all; 
    h(2) = loglog(labels(i),anN_2(i),'ko','MarkerSize',12,'LineWidth',3); 
    hold all;
    h(3)=loglog(labels(i),anN_3(i),'kx','MarkerSize',14,'LineWidth',3); 
    hold all;  
end
%Finalize k figure
figure(k_fig); 
%Plot bounds from previous
hold on; 
semilogx([1 1e7],[26 26],'--','LineWidth',2,'Color',[0.5 0.5 0.5]); 
hold on; 
semilogx([1 1e7], [36 36],'--','LineWidth',2,'Color',[0.5 0.5 0.5]); 
xlabel('Degradation time [s]','FontSize',30); 
ylabel('k [-]','FontSize',30); 
axis([0 1e7 0 50]);
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
%Finalize tau_n0 figure
figure(tau_n0_fig); 
legend(h,'Ti','Mo','W');
set(get(h(1),'Parent'),'YScale','log');
xlabel('Degradation time [s]','FontSize',30); 
ylabel('N_t [cm^{-3}]','FontSize',30); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1e7 1e10 1e13]);

%% Look at how tau_n0 changes with sensitivity
clear all; close all; 
%Manually enter fit values
best_fit = {[185.887526819101;5.14136156477898],...%original
            [167.862890697482; 5.16203332822537],...%increase injection
            [213.894670034334; 5.00648345496999],...%decrease injection
            [203.828608040808; 5.69833765958448],... %increase lifetime
            [172.121993432101; 4.50217342989029],... %decrease lifetime
            [192.937816540636; 5.09677150266217],...%increase Joe
            [200.242867847582; 4.77821069627767]} %decrease Joe
label = {'original','increase inj','decrease inj','increase tau','decrease tau','increase Joe','decrease Joe'}; 
doping = 9.09e15; 
T = 300; 
type = 'p'; 
%hard code the parameters as used for publication
n0 = 1.02e4; 
p0 =  9.09e15; 
NC = 3e19; 
NV = 1e19; 
vth_e = 2.05e7; 
vth_h = 1.69e7; 
Eg = 1.1242; 
k_B = 8.61733238e-5;
defect1 = figure;
tau_defect1 = figure; 
[m,n] = size(best_fit); 
Et= cell(m,n); 
k= cell(m,n); 
alphanN= cell(m,n); 
for i = 1:n
    best_fit_hold = best_fit{i};
    [fits,num_defect] = size(best_fit_hold); 
    for j = 1:num_defect
        best_fit_now = best_fit_hold(:,j); 
        %Define the energy levels for evaluation
        Et_now = linspace(0,Eg,250); %eV
        Q = zeros(size(Et_now)); 
        alphanN_now = zeros(size(Et_now));
        k_now = zeros(size(Et_now)); 
        A = best_fit_now(1)+best_fit_now(2); %X -> 1
        B = best_fit_now(2); %X -> 0
        C = best_fit_now(1)/A; %slope/X -> 1
        for l = 1:length(Et_now)
            %Calculate n1
            n1 = NC*exp(-(Eg-Et_now(l))/(k_B*T)); 
            %Calculate p1
            p1 = NV*exp(-Et_now(l)/(k_B*T)); 
            %Calculate the Q values for these defects
            Q(l) = (C+(p1/p0))/(1-(n1/p0)-C);
            %Calculate the quantity alphan*Nt for these defects
            alphanN_now(l) = (1/B)*(1+((1/p0)*((Q(l)*n1)+p1))); 
            %Calculate the k values for these defects
            k_now(l) = Q(l)*vth_h/vth_e; 
        end
        %Get rid of any negative k values
%         indices = find(k_now<0); 
%         Et_now(indices) = []; 
%         k_now(indices) = []; 
%         alphanN_now(indices) = []; 
        Et{i,j} = Et_now; 
        k{i,j} = k_now; 
        alphanN{i,j} = alphanN_now; 
    end
    figure(defect1); 
    h1(i)=plot(Et{i,1},k{i,1},'-','LineWidth',3); 
%     plot(Et{i,1},k{i,1},'-','LineWidth',3,'Color',co{i}); 
%     label(i,1) = T; 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',3);
%     plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',3,'Color',co{i});
    hold all;
end
figure(defect1); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('k [-]','FontSize',30);
legend(h1,label');
title('Defect 1','FontSize',30); 
figure(tau_defect1); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('\tau_{n0} [\mus]','FontSize',30);
legend(h3,label');
title('Defect 1','FontSize',30);

%Let's calculate the error and save it
original = 1; 
min = 5;%decrease tau
max = 4;%increase tau
alphanN_up = (1./alphanN{max,1})-(1./alphanN{original,1}); 
alphanN_down = (1./alphanN{original,1})-(1./alphanN{min,1}); 
alphanN_up = alphanN_up./(1./alphanN{original,1}); 
alphanN_down = alphanN_down./(1./alphanN{original,1}); 
Et = Et{original,1}; 
% save('Estimated_taun0error_from_JPV.mat','alphanN_up','alphanN_down','Et');

%% Make values for plotting in origin: DEGRADATION
clear all; close all; clc;
%Load the errors and save 
error_k = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\Publication\Estimated_kerror_from_JPV.mat';
load(error_k); 
Et_errork = Et; 
k_up = up; 
k_down = down; 
error_taun0 = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\Publication\Estimated_taun0error_from_JPV.mat';
load(error_taun0); 
Et_errortau = Et; 
tau_up = alphanN_up; 
tau_down = alphanN_down; 
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\Processed degraded';
load([directory '\fitted_defect_parameters_deg.mat']); 
labels = [500 1000 10000 100000 200021 300297 380757]; %degradation
%Now plot the k-values and tau_n0 values for each defect at midgap
[num_meas,num_defects] = size(k); 
tau_store = zeros(num_meas,3); %for Tidd, Mo, W
k_store = zeros(num_meas,3);
k_upper = zeros(num_meas,3);
k_lower = zeros(num_meas,3);
tau_upper = zeros(num_meas,3);
tau_lower = zeros(num_meas,3);
k_mg = zeros(num_meas,1); 
k_mg_upper = zeros(num_meas,1); 
k_mg_lower = zeros(num_meas,1); 
for i = 2 %we only care about defect 2
    midgap = 1.124/2; 
    tau_n0_fig = figure;
    k_fig = figure;
    for j = 1:num_meas
        Et_now = Et{j,i};
        mg_index = find(min(abs(Et_now-midgap))==abs(Et_now-midgap)); %closest index to midgap
        k_now = k{j,i}; 
        alphanN_now = alphanN{j,i};
        k_mg(j,1) = k_now(mg_index);
        k_mg_upper(j,1) = k_up(find(min(abs(Et_errork-midgap))==abs(Et_errork-midgap)));
        k_mg_lower(j,1) = k_down(find(min(abs(Et_errork-midgap))==abs(Et_errork-midgap)));
        k_store(j,1) = k_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        k_store(j,2) = k_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        k_store(j,3) = k_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
        tau_store(j,1) = 1/alphanN_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        tau_store(j,2) = 1/alphanN_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        tau_store(j,3) = 1/alphanN_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
        kerror_upper(j,1) = k_up(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        kerror_upper(j,2) = k_up(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        kerror_upper(j,3) = k_up(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W
        kerror_lower(j,1) = k_down(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        kerror_lower(j,2) = k_down(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        kerror_lower(j,3) = k_down(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W
        tauerror_upper(j,1) = tau_up(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        tauerror_upper(j,2) = tau_up(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        tauerror_upper(j,3) = tau_up(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W
        tauerror_lower(j,1) = tau_down(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        tauerror_lower(j,2) = tau_down(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        tauerror_lower(j,3) = tau_down(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W  
    end
end
%Let's try plotting for different possible defects
sigman_1 = 1.5e-15; %Tidd
sigman_2 = 1.6e-14; %Mo
sigman_3 = 1.7e-14; %W
%Modify the units so that they are in seconds rather than microseconds
alphanN_store = (1./tau_store).*1e6; 
alphanN_upper = (1./(tau_store+(tauerror_upper.*tau_store))).*1e6; 
alphanN_lower = (1./(tau_store-(tauerror_lower.*tau_store))).*1e6; 
vth_e = 2.05e7; %cm/s
alphanN_store = alphanN_store./vth_e; 
alphanN_upper = alphanN_upper./vth_e; 
alphanN_lower = alphanN_lower./vth_e; 
alphanN_store(:,1) = alphanN_store(:,1)./sigman_1; 
alphanN_upper(:,1) = alphanN_upper(:,1)./sigman_1; 
alphanN_lower(:,1) = alphanN_lower(:,1)./sigman_1; 
alphanN_store(:,2) = alphanN_store(:,2)./sigman_2; 
alphanN_upper(:,2) = alphanN_upper(:,2)./sigman_2; 
alphanN_lower(:,2) = alphanN_lower(:,2)./sigman_2; 
alphanN_store(:,3) = alphanN_store(:,3)./sigman_3; 
alphanN_upper(:,3) = alphanN_upper(:,3)./sigman_3; 
alphanN_lower(:,3) = alphanN_lower(:,3)./sigman_3; 
alphanN_upper = alphanN_store-alphanN_upper;
alphanN_lower = alphanN_lower-alphanN_store; 
k_upper = (k_store.*kerror_upper); 
k_lower = (k_store.*kerror_lower); 
k_mg_upper = (k_mg.*k_mg_upper); 
k_mg_lower = (k_mg.*k_mg_lower); 

%% Make values for plotting in origin: REGENERATION
clear all; close all; clc;
%Load the errors and save 
error_k = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\Publication\Estimated_kerror_from_JPV.mat';
load(error_k); 
Et_errork = Et; 
k_up = up; 
k_down = down; 
error_taun0 = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\Publication\Estimated_taun0error_from_JPV.mat';
load(error_taun0); 
Et_errortau = Et; 
tau_up = alphanN_up; 
tau_down = alphanN_down; 
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\Processed regeneration\Excluding low injection data';
load([directory '\fitted_defect_parameters_regen.mat']); 
labels = [380757 826037 1105807 1530090 2061150 2325570 2574840 2733780]; %regeneration
%Now plot the k-values and tau_n0 values for each defect at midgap
[num_meas,num_defects] = size(k); 
tau_store = zeros(num_meas,3); %for Tidd, Mo, W
k_store = zeros(num_meas,3);
k_upper = zeros(num_meas,3);
k_lower = zeros(num_meas,3);
tau_upper = zeros(num_meas,3);
tau_lower = zeros(num_meas,3);
k_mg = zeros(num_meas,1); 
k_mg_upper = zeros(num_meas,1); 
k_mg_lower = zeros(num_meas,1); 
for i = 2 %we only care about defect 2
    midgap = 1.124/2; 
    tau_n0_fig = figure;
    k_fig = figure;
    for j = 1:num_meas
        Et_now = Et{j,i};
        mg_index = find(min(abs(Et_now-midgap))==abs(Et_now-midgap)); %closest index to midgap
        k_now = k{j,i}; 
        alphanN_now = alphanN{j,i};
        k_mg(j,1) = k_now(mg_index);
        k_mg_upper(j,1) = k_up(find(min(abs(Et_errork-midgap))==abs(Et_errork-midgap)));
        k_mg_lower(j,1) = k_down(find(min(abs(Et_errork-midgap))==abs(Et_errork-midgap)));
        k_store(j,1) = k_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        k_store(j,2) = k_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        k_store(j,3) = k_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
        tau_store(j,1) = 1/alphanN_now(find(min(abs(Et_now-0.259))==abs(Et_now-0.259))); %Tidd
        tau_store(j,2) = 1/alphanN_now(find(min(abs(Et_now-0.28))==abs(Et_now-0.28))); %Mo
        tau_store(j,3) = 1/alphanN_now(find(min(abs(Et_now-0.55))==abs(Et_now-0.55))); %W
        kerror_upper(j,1) = k_up(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        kerror_upper(j,2) = k_up(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        kerror_upper(j,3) = k_up(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W
        kerror_lower(j,1) = k_down(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        kerror_lower(j,2) = k_down(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        kerror_lower(j,3) = k_down(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W
        tauerror_upper(j,1) = tau_up(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        tauerror_upper(j,2) = tau_up(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        tauerror_upper(j,3) = tau_up(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W
        tauerror_lower(j,1) = tau_down(find(min(abs(Et_errork-0.259))==abs(Et_errork-0.259))); %Tidd
        tauerror_lower(j,2) = tau_down(find(min(abs(Et_errork-0.28))==abs(Et_errork-0.28))); %Mo
        tauerror_lower(j,3) = tau_down(find(min(abs(Et_errork-0.55))==abs(Et_errork-0.55))); %W  
    end
end
%Let's try plotting for different possible defects
sigman_1 = 1.5e-15; %Tidd
sigman_2 = 1.6e-14; %Mo
sigman_3 = 1.7e-14; %W
%Modify the units so that they are in seconds rather than microseconds
alphanN_store = (1./tau_store).*1e6; 
alphanN_upper = (1./(tau_store+(tauerror_upper.*tau_store))).*1e6; 
alphanN_lower = (1./(tau_store-(tauerror_lower.*tau_store))).*1e6;  
vth_e = 2.05e7; %cm/s
alphanN_store = alphanN_store./vth_e; 
alphanN_upper = alphanN_upper./vth_e; 
alphanN_lower = alphanN_lower./vth_e; 
alphanN_store(:,1) = alphanN_store(:,1)./sigman_1; 
alphanN_upper(:,1) = alphanN_upper(:,1)./sigman_1; 
alphanN_lower(:,1) = alphanN_lower(:,1)./sigman_1; 
alphanN_store(:,2) = alphanN_store(:,2)./sigman_2; 
alphanN_upper(:,2) = alphanN_upper(:,2)./sigman_2; 
alphanN_lower(:,2) = alphanN_lower(:,2)./sigman_2; 
alphanN_store(:,3) = alphanN_store(:,3)./sigman_3; 
alphanN_upper(:,3) = alphanN_upper(:,3)./sigman_3; 
alphanN_lower(:,3) = alphanN_lower(:,3)./sigman_3; 
alphanN_upper = alphanN_store-alphanN_upper;
alphanN_lower = alphanN_lower-alphanN_store;
k_upper =(k_store.*kerror_upper); 
k_lower = (k_store.*kerror_lower); 
k_mg_upper = (k_mg.*k_mg_upper); 
k_mg_lower = (k_mg.*k_mg_lower); 


