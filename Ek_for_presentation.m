%Ek plots for presentation 
clear all; close all; clc;
dirname = 'C:\Users\Mallory Jensen\Dropbox (MIT)\TIDLS data - UROP\Data\analysis by ingot\22\22 top analysis'; 
load([dirname '\E-k.mat']); 
load([dirname '\lifetime_breakdown.mat']); 
load([dirname '\meas_info.mat']);
% indices = [1 2 3]; 
% labels = {'-75C';'-25C';'25C'}; 
% indices = [1 2]; 
% labels = {'25C','100C'};
% indices = [2 3 4]; 
% labels = {'-25C','25C','100C'};
% indices = [ 1 2 3 4];
% labels = {'-75C';'-25C';'25C';'100C'};
indices = [2 3]; 
labels = {'-25C','25C'};


[Efi,Efv,p0,n0,Eiv] = adv_Model_gen(25+273.15,info(1).doping,'n');
defect_file = 'C:\Users\Mallory Jensen\Dropbox (MIT)\TIDLS data - UROP\Data\analysis by ingot\defects.xlsx';
[num,txt]=xlsread(defect_file); 
names = txt(2:end,1);
Etv = num(:,1); 
kval = num(:,3); 
Eti = Etv-Eiv; 

%colors
coldest = [0 0 1];
cold = [0 0.75 0.75];
room = [0.93 0.69 0.13];
warm = [0.83 0.33 0.1];
% colors = {coldest,cold,room}; 
% colors = {room, warm};
% colors = {coldest,cold, room, warm};
colors = {cold,room};

%graph properties
lw = 3; 
lwgraph = 2; 
ticklengths = [0.025; 0.05]; 
xlimits = [min(Et_vector) max(Et_vector)];
xlim_lifetime = [1e13 1e16];

SRHtau = figure;
Ek = figure;
Etau = figure;
stdev = figure; 
avg = figure;
stdev_tau = figure;
avg_tau = figure; 

std_dev = zeros(100,1); 
average = zeros(100,1); 
std_dev(:,1) = NaN; average(:,1) = NaN;
std_dev_tau = zeros(100,1); average_tau = zeros(100,1); 
std_dev_tau(:,1) = NaN; average_tau(:,1) = NaN;

for i = 1:length(indices)
    figure(SRHtau); 
    loglog(lifetime_breakdown(indices(i)).deltan,...
        lifetime_breakdown(indices(i)).tau_SRH.*1e6,'LineWidth',lw,'Color',...
        colors{i}); 
    hold all; 
    
    figure(Ek); 
    Et_now = Et{indices(i),1}; k_now = k{indices(i),1}; 
    semilogy(Et_now,k_now,'LineWidth',lw,'Color',colors{i});
    hold all; 
    
    figure(Etau); 
    tau_now = alphanN{indices(i),1}; tau_now = (1./tau_now).*1e6; 
    semilogy(Et_now,tau_now,'LineWidth',lw,'Color',colors{i});
    hold all;  
end
    
for i = 1:length(Et_vector)
    for j = 1
        for q = 1:length(indices)
            k_now = k{indices(q),j};
            tau_now = alphanN{indices(q),j}; tau_now = (1./tau_now).*1e6; 
            if isnan(k_now(i))==1
               all_T(j,q) = NaN; 
            else
                all_T(j,q) = k_now(i);
            end
            if isnan(tau_now(i))==1
               all_T_tau(j,q) = NaN; 
            else
                all_T_tau(j,q) = tau_now(i);
            end
        end
        std_dev(i,j) = nanstd(all_T(j,:)); 
        %We should also track the mean
        average(i,j) = nanmean(all_T(j,:));
        std_dev_tau(i,j) = nanstd(all_T_tau(j,:)); 
        average_tau(i,j) = nanmean(all_T_tau(j,:)); 
    end
end

figure(stdev);
semilogy(Et_vector,std_dev,'k','LineWidth',lw);
set(gca,'FontSize',18); 
xlabel('E_t-E_i [eV]','FontSize',20);
ylabel('standard deviation in \it{k}','FontSize',20); 
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlimits);
hgsave(stdev,[dirname '\Et-kstdev']); 

figure(avg); 
semilogy(Et_vector,average,'k','LineWidth',lw); 
set(gca,'FontSize',18); 
xlabel('E_t-E_i [eV]','FontSize',20);
ylabel('average \it{k}','FontSize',20);
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlimits);
hgsave(avg,[dirname '\Et-kavg']);

figure(stdev_tau);
semilogy(Et_vector,std_dev_tau,'k','LineWidth',lw);
set(gca,'FontSize',18); 
xlabel('E_t-E_i [eV]','FontSize',20);
ylabel('standard deviation in \tau_p_0','FontSize',20); 
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlimits);
hgsave(stdev_tau,[dirname '\Et-taustdev']); 

figure(avg_tau); 
semilogy(Et_vector,average_tau,'k','LineWidth',lw); 
set(gca,'FontSize',18); 
xlabel('E_t-E_i [eV]','FontSize',20);
ylabel('average \tau_p_0','FontSize',20);
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlimits);
hgsave(avg_tau,[dirname '\Et-tauavg']);

figure(SRHtau); 
set(gca,'FontSize',18); 
xlabel('excess carrier density [cm^-^3]','FontSize',20); 
ylabel('lifetime [\mus]','FontSize',20);
legend(labels); 
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlim_lifetime);
hgsave(SRHtau,[dirname '\SRH_lifetimes']);

figure(Ek); 
set(gca,'FontSize',18); 
xlabel('E_t-E_i [eV]','FontSize',20);
ylabel('{\itk} [-]','FontSize',20); 
legend(labels); 
hold all;
semilogy(Eti,kval,'r.','MarkerSize',20); 
text(Eti+0.01,kval,names);
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlimits);
hgsave(Ek,[dirname '\Et-k']);

figure(Etau);
set(gca,'FontSize',18); 
xlabel('E_t-E_i [eV]','FontSize',20);
ylabel('\tau_p_0 [\mus]','FontSize',20); 
legend(labels); 
ax = gca; 
ax.LineWidth = lwgraph; 
ax.TickLength = ticklengths;
xlim(xlimits);
hgsave(Etau,[dirname '\Et-tau']);