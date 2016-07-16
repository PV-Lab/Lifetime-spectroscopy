clear all; close all; 

sample_nos = {'19-2','19-3','19-5_1-1','19-5_1-64','19-6',...
    '20-1','20-2','20-4','20-6',...
    '123-2_1-1','123-2_1-64','123-3_1-1','123-3_1-64','123-5_1-1','123-5_1-64','123-6_1-1','123-6_1-64',...
    '124-1_1-1','124-1_1-64','124-2','124-4_1-1','124-4_1-64','124-6_1-1','124-6_1-64',...
    '16-2-18-N'};

for i = 1:length(sample_nos)
    sample_no = sample_nos{i};
    
    %Get the data from Excel
    dirname = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 3\4 days after anneal\Retaking averages\' sample_no];
    saveName = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 3\4 days after anneal\Retaking averages\' sample_no '\all_XLS_data.mat'];
    process_xls_data(dirname,saveName);
    
    filename = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 3\4 days after anneal\Retaking averages\' sample_no '\all_XLS_data.mat'];

    saveStart = ['C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 3\4 days after anneal\Retaking averages\' sample_no '\'];
    [tau_mean,deltanq]=average_tau(filename,saveStart);
    saveFile = [saveStart sample_no '_averageTau.mat'];
    save(saveFile,'tau_mean','deltanq');
    
end