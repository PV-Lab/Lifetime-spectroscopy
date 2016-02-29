clear all; close all; 

sample_nos = {'19-1','19-2','19-3','19-4','19-5','19-6'...
    '20-1','20-2','20-3','20-4','20-5','20-6',...
    '123-1 1-64','123-2','123-3','123-4 1-1','123-5','123-6',...
    '124-1','124-2','124-3 1-64','124-4','124-5','124-6'};
N_dop = [2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,2.8e15,4.1e15,4.1e15,4.1e15,4.1e15,4.1e15,4.1e15,3.8e15,3.8e15,3.8e15,3.8e15,3.8e15,3.8e15];
type = 'n';
T = 300; 
W = [170,180,180,170,180,180,180,180,170,180,170,180,170,180,180,170,180,180,180,180,170,180,170,180]; %cm
D = [11.84,11.84,11.84,11.84,11.84,11.84,11.84,11.84,11.84,11.84,11.84,11.84,11.72,11.72,11.72,11.72,11.72,11.72,11.75,11.75,11.75,11.75,11.75,11.75]; %cm2/s
SRV_file = 'C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 1\No box\15-9-21-N\SRV_deltan.mat';
deltanq = logspace(14,17,500);

for i = 1:length(sample_nos)
    D_now = D(i);
    W_now = W(i); 
    N_dop_now = N_dop(i);
    sample_no = sample_nos{i};
    
    measured_file = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 1\No box\' sample_no '\LifetimeData.mat'];
    saveFile = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 1\No box\' sample_no '\' sample_no '_tauBreakdown.mat'];
    saveStart = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 1\No box\' sample_no '\' sample_no];
    [tau_SRH,deltanq] = calculate_SRH(measured_file,SRV_file,deltanq,N_dop_now,T,type,W_now,D_now,saveFile,saveStart,sample_no);
end