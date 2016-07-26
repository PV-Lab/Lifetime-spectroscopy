clear all; close all; clc; 
%Run through different temperatures and evaluate optical constant
dirname = 'C:\Users\Mallory\Documents\Australia\Optical constant\';
%folder names
folders = {'25C','75C','100C','175C','200C'}; 
T = [25 75 100 175 200]; 
OC_range = [0.8 1.2]; 
resolution = 20; 
inj_range = [3e15 1e16]; %assume this can be the same for all temperatures
for i = 1:length(folders)
    location = [dirname folders{i}];
    [fileList,fileListShort] = getAllFiles(location);
    %gen is first (1-1)
    filename_gen = fileList{1}; 
    filename_trans = fileList{2}; 
    [deltan_trans,tau_trans,dataSave_gen,OC_values,MSE,min_index]=modify_OC(filename_trans,filename_gen,OC_range,inj_range,resolution);
    OC_win(i) = OC_values(min_index); 
end

figure; 
plot(T,OC_win,'o'); 
xlabel('Temperature [C]','FontSize',20); 
ylabel('Best optical constant [-]','FontSize',20); 
title('Float zone test sample','FontSize',20); 
