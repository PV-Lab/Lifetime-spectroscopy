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
measurement_log = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\from SERIS\Symmetrical_Asymmetrical Lifetime samples\measurement_summary.xlsx';
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\from SERIS\Symmetrical_Asymmetrical Lifetime samples\Symmetrical-AsymmetricalLifetime';
%Define the samples as they are listed in the filenames
samples = {'Assymetric','Symmetric'};
%Read in the data with the times
filename_details = cell(length(samples),2); 
for i = 1:length(samples)
    [num,txt,raw] = xlsread(measurement_log,samples{i}); 
    filename_details{i,1} = txt(2:end,1:2); 
    filename_details{i,2} = num; 
end
filename_start = [directory '\']; 
filename_end = {'.xlsm' '.xlsm'}; 

%Time intervals
% times = [0, 10:10:100, 200:100:1000, 2000:1000:10000 20000:10000:50000];

colors = {'r','g','b','m','c','y'};

zero_filenames = {'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\from SERIS\Symmetrical_Asymmetrical Lifetime samples\Symmetrical-AsymmetricalLifetime\Assymetric_0hrs.xlsm',...
    'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\from SERIS\Symmetrical_Asymmetrical Lifetime samples\Symmetrical-AsymmetricalLifetime\Symmetric_0hrs.xlsm'};

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
            filename = [filename_start  samples{i} '_' times_str{j} filename_end{i}];
        end
        %For the WCT-120TS spreadsheet
        data = xlsread(filename,'RawData','G4:I124');
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
        h(j)=loglog(deltan,tau.*1e6,'-','LineWidth',3);
        label(j) = times_num(j); 
        hold all; 
%         if round(log10(times_num(j)))==log10(times_num(j))
%             h(count) = loglog(deltan,tau.*1e6,[colors{count} '--'],'LineWidth',3);
%             label(count) = times_num(j); 
%             count = count+1;  
%         else
%             loglog(deltan,tau.*1e6,'.','MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',6);
%         end
%         hold on; 
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
injection = 8e14; %cm-3
% injection = 1e15; 
figure; 
for i = 1:length(samples)
    %Get the number of files
    times_str = filename_details{i,1}; 
    times_str = times_str(:,2); 
    times_num = filename_details{i,2}; 
    times_num = times_num(:,1); 
    for j = 1:length(times_num)
        data_now = dataSave{i,j}; 
        lifetime_deg(i,j) = interp1(data_now(:,1),data_now(:,2),injection); 
    end
    semilogx(times_num,lifetime_deg(i,1:length(times_num)).*1e6,'o','MarkerSize',12,'LineWidth',2); 
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
    times_num = times_num(:,1); 
%     lifetime_deg_norm(i,:) = lifetime_deg(i,:)./max(lifetime_deg(i,:)); 
    lifetime_deg_norm(i,:) = lifetime_deg(i,:)./lifetime_deg(i,1);
    semilogx(times_num,lifetime_deg_norm(i,1:length(times_num)),'o','MarkerSize',12,'LineWidth',2); 
    hold all;
end
xlabel('Degradation time [s]','FontSize',30); 
ylabel('Normalized lifetime [-]','FontSize',30);
legend(samples); 

%Save the data 
save([directory '\processed_data_8e14_20170222.mat'],'dataSave','lifetime_deg_norm','lifetime_deg','filename_details','samples','times_num');

%% Assess lifetime fitting with same parameters as previous publication
clc;clear all; close all; 
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\from SERIS\Symmetrical_Asymmetrical Lifetime samples\Symmetrical-AsymmetricalLifetime';
load([directory '\processed_data_8e14_20170222.mat']);

sample_index = 1; %asymmetric
to_calc = [2:1:13]; %asymmetric
% sample_index = 2; %symmetric
% to_calc = [2:1:12]; %symmetric

%Find the fully degraded state which will correspond to the minimum
%normalized degraded lifetime
max_deg_index = find(lifetime_deg_norm(sample_index,:)==min(lifetime_deg_norm(sample_index,1:12))); 
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
[num_samples,num_measurements] = size(dataSave); 
times_num = filename_details{sample_index,2}; times_num = times_num(:,4); 
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
title('Asymmetric','FontSize',20); 
legend(h,num2str(labels')); 
print(SRHfig,'-dpng','-r0',[directory '\SRH_lifetimes_asymmetric.png']);
hgsave(SRHfig,[directory '\SRH_lifetimes_asymmetric']);
figure(taufig); 
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('Lifetime [s]','FontSize',20);
title('Asymmetric','Fontsize',20);
legend(h2,num2str([0;labels'])); 
print(taufig,'-dpng','-r0',[directory '\Measured_lifetimes_asymmetric.png']);
hgsave(taufig,[directory '\Measured_lifetimes_asymmetric']);
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
    xlswrite([directory '\Linearized_data_symmetric_3.xlsx'],[X,tau_rev],['Sheet' num2str(i)]); 
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
    xlswrite([directory '\Linearized_data_symmetric_3.xlsx'],to_write,['Sheet' num2str(i)],'C1:E6'); 
end
%Now pause and refine fits in Excel. 
%% After fits refined in Excel, make Ek curves
clear all; close all; 
directory = 'C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\from SERIS\Symmetrical_Asymmetrical Lifetime samples\Symmetrical-AsymmetricalLifetime';
% to_calc = [2:1:13]; %asymmetric
% label = {'61200','136800','219600','298800','381600','399600','471600','554400','676800','835200','842400','892800'};
to_calc = [2,4,5,7,8,9,10]; %symmetric, excluding weird fits
label = {'136800','298800','381600','471600','554400','576000','676800'};
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
filename = [directory '\IDLS two and three curve fitting - symmetric.xlsm'];
% fits = xlsread(filename,'Summary','B3:E14'); %They should be in the same order as those in the MATLAB file
fits = xlsread(filename,'for Ek','B3:E9');
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
    h1(i)= plot(Et{i,1},k{i,1},'-','LineWidth',3); 
    hold all; 
    figure(tau_defect1); 
    h3(i)=plot(Et{i,1},1./alphanN{i,1},'-','LineWidth',3);
    hold all;
    figure(defect2);
    h2(i)=plot(Et{i,2},k{i,2},'-','LineWidth',3);
    hold all;
    figure(tau_defect2); 
    h4(i)=plot(Et{i,2},1./alphanN{i,2},'-','LineWidth',3); 
    hold all; 
end
figure(defect1); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('k [-]','FontSize',30);
legend(h1,label);
%Let's also plot the identified range from previous paper 26-36 (k)
hold all; 
plot([0 1.13],[26 26],'--','LineWidth',2,'Color',[0.5 0.5 0.5]); 
hold all; 
plot([0 1.13],[36 36],'k--','LineWidth',2,'Color',[0.5 0.5 0.5]);
title('Defect 1','FontSize',30); 
figure(defect2);
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 50]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('k [-]','FontSize',30);
legend(h2,label);
title('Defect 2','FontSize',30); 
figure(tau_defect1); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('\tau_{n0} [\mus]','FontSize',30);
legend(h3,label);
title('Defect 1','FontSize',30);
figure(tau_defect2); 
set(gca,'FontSize',20);
set(gca,'LineWidth',2);
axis([0 1.124 0 100]);
xlabel('E_t-E_v [eV]','FontSize',30); 
ylabel('\tau_{n0} [\mus]','FontSize',30);
legend(h4,label);
title('Defect 2','FontSize',30);
save([directory '\fitted_defect_parameters_symmetric.mat'],'Et','alphanN','k','two_defects','T','doping','type','label');

