%Plot the degradation curves of all samples
clc;clear all; close all; 
measurement_log = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\measurement_summary.xlsx';
directory = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC';
%Define the samples as they are listed in the filenames
samples = {'64-5' '69-5','FZ'};
%Read in the data with the times
filename_details = cell(length(samples),2); 
for i = 1:length(samples)
    [num,txt,raw] = xlsread(measurement_log,samples{i}); 
    filename_details{i,1} = txt(2:end,1:2); 
    filename_details{i,2} = num; 
end
filename_start = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\'; 
filename_end = {'_1-64_avg5.xlsm' '.xlsm' '_avg5.xlsm'}; 

%Time intervals
times = [0, 10:10:100, 200:100:1000, 2000:1000:10000 20000:10000:50000];

colors = {'r','g','b','m','c','y'};

zero_filenames = {'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\64-5\64-5_beforeDeg_1-64_avg5.xlsm' 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\69-5\69-5_afterAnneal_avg5.xlsm' 'C:\Users\Mallory\Documents\PERC mc-Si degradation\Experiment 0\rev resistivity, OC\FZ\FZ_afterAnneal_avg5.xlsm'}; 

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
            %Make the filename the way we expect it
            filename = [filename_start samples{i}  '\' samples{i} '_' times_str{j} filename_end{i}];
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

%Choose one injection level and plot the lifetime at that injection level
%Now that we have data stored, it's easy enough just to cycle through it
injection = 5e14; %cm-3
figure; 
for i = 1:length(samples)
    %Get the number of files
    times_str = filename_details{i,1}; 
    times_str = times_str(:,2); 
    times_num = filename_details{i,2}; 
    times_num = times_num(:,4); 
    for j = 1:length(times_num)
        data_now = dataSave{i,j}; 
        lifetime_deg(i,j) = interp1(data_now(:,1),data_now(:,2),injection); 
    end
    plot(times_num,lifetime_deg(i,1:length(times_num)).*1e6,'o','MarkerSize',12,'LineWidth',2); 
    hold all; 
end
xlabel('Degradation time [s]','FontSize',30); 
ylabel('Lifetime [\mus]','FontSize',30);
legend(samples); 

%Normalized degradation
figure; 
for i = 1:length(samples); 
    %Get the number of files
    times_str = filename_details{i,1}; 
    times_str = times_str(:,2); 
    times_num = filename_details{i,2}; 
    times_num = times_num(:,4); 
    lifetime_deg_norm(i,:) = lifetime_deg(i,:)./max(lifetime_deg(i,:)); 
    plot(times_num,lifetime_deg_norm(i,1:length(times_num)),'o','MarkerSize',12,'LineWidth',2); 
    hold all;
end
xlabel('Degradation time [s]','FontSize',30); 
ylabel('Normalized lifetime [-]','FontSize',30);
legend(samples); 

%Save the data 
save([directory '\processed_data_20160811.mat'],'dataSave','lifetime_deg_norm','filename_details','samples');

%% Given the loaded data, try now to analyze the evolution along the degradation curve
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
save([directory '\fitted_defect_parameters.mat'],'Et','alphanN','k','two_defects','T','doping','type','samples','labels');
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
        k_store(j,i) = k_now(mg_index);
        alphanN_store(j,i) = alphanN_now(mg_index); 
    end
end
%I only care about defect 2 for the moment
k_def2 = k_store(:,2); 
alphanN_def2 = alphanN_store(:,2); 
%For k, we define a general error based on previous experiments
k_error = k_def2.*0.1; %+/- 10%
figure(k_fig); 
colors = {[0 0.4980 0],[1 0 0],[0 0.749 0.749],[0.749 0 0.749],[0.749 0.749 0],[0.2471 0.2471 0.2471]};
for i = 1:length(k_def2)
    errorbar(labels(i)',k_def2(i),k_error(i),'x','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
    hold on;
end
xlabel('Degradation time [s]','FontSize',30); 
ylabel('\sigma_n/\sigma_p [-]','FontSize',30); 
axis([0 1e6 0 50]);
%For tau_n0, we can define the error based on perturbations in the measured
%lifetime
error_upper = [1.13 1.09 1.07 1.06 1.06 1.06];
error_lower = [0.82 0.9 0.93 0.93 0.93 0.93]; 
%Let's try plotting for different possible defects
sigman_1 = 1e-13; 
sigman_2 = 1e-14; 
sigman_3 = 1e-15; 
%Modify the units so that they are in seconds rather than microseconds
alphanN_def2 = alphanN_def2.*1e6; %1/s
%Let's take out the thermal velocity for clarity. These terms are
%Nt*sigma_n*vth_e
vth_e = 2.05e7; %cm/s
alphanN_def2 = alphanN_def2./vth_e; %Nt*sigma_n, units cm-1
%Now also take out the different sigma values
anN_1 = alphanN_def2./sigman_1; 
anN_2 = alphanN_def2./sigman_2; 
anN_3 = alphanN_def2./sigman_3; 
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
    h(1)=loglog(labels(i),anN_1(i),'x','MarkerSize',12,'LineWidth',2,'Color',colors{i}); 
    hold all;
    h(2)=loglog(labels(i),anN_3(i),'s','MarkerSize',12,'LineWidth',2,'Color',colors{i}); 
%     h(1)=errorbar(labels(i),anN_1(i),L_1(i),U_1(i),'x','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
%     hold on;
%     errorbar(labels(i),anN_2(i),L_2(i),U_2(i),'o','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
%     hold all;
%     h(2)=errorbar(labels(i),anN_3(i),L_3(i),U_3(i),'s','MarkerSize',12,'LineWidth',3,'Color',colors{i}); 
end
legend(h,'10^{-13}','10^{-15}');
set(get(h(1),'Parent'),'YScale','log');
xlabel('Degradation time [s]','FontSize',30); 
ylabel('N_t [cm^{-3}]','FontSize',30); 
        
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


    