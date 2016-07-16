%Comparing SRH lifetimes
saveStart = 'C:\Users\Mallory\Documents\Sinton visit\May 1 2016\124-1';
filename1 = 'C:\Users\Mallory\Documents\Sinton visit\May 1 2016\124-1\lifetime_breakdown.mat'; 
filename2 = 'C:\Users\Mallory\Documents\Non-contact crucible\Australia experiments\June 22 2016\124-1\lifetime_breakdown.mat'; 

load(filename1); 

figure; 
loglog(deltanq,tau_SRH.*1e6);
tau_SRH1 = tau_SRH; 
load(filename2); 
hold all; 
loglog(deltanq,tau_SRH.*1e6);
tau_SRH2 = tau_SRH; 
legend('Measurement 1','Measurement 2'); 

xlabel('Excess carrier density (cm^-^3)','FontSize',30);
ylabel('Lifetime (\mus)','FontSize',30);
axis([5e13 1e17 0 10000]);
% print('-dpng','-r0',[saveStart '\SRH_comparison.png']);

figure;
semilogx(deltanq,((tau_SRH1-tau_SRH2)./tau_SRH1).*100); 