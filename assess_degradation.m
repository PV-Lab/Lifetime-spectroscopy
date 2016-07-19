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
save([directory '\processed_data.mat'],'dataSave','lifetime_deg_norm','filename_details','samples');

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
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h3,num2str(label));
title('Defect 1','FontSize',30);
figure(tau_defect2); 
xlabel('E_t-E_v [eV]','FontSize',20); 
ylabel('\tau_{n0} [s]','FontSize',20);
legend(h4,num2str(label));
title('Defect 2','FontSize',30);
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
times = [10000 20000 100000 200000 470000 610000 870000]; 
folders = {'10000s' '20000s' '100000s' '200000s' '470000s' '610000s' '870000s'}; 
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
        filename = [filename_start '\' folders{j} '\' samples{i} '_' num2str(times(j)) '_' num2str(exposure(i)) 's_' num2str(LP) filename_end]; 
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


    