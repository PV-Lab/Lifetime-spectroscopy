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
%% First process the raw data
clear all; close all; clc; 
dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 17 2017'; 
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','C-1','C-2','H-1','H-2','FZ'};
samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12'};
for index = 1:length(samples)
    [fileList,fileListShort] = getAllFiles([dirname '\' samples{index}]); 
    savename = [dirname '\' samples{index} '\Raw_data.mat']';
    process_xls_data([dirname '\' samples{index}],savename);
    %We also want to store all of the information for each file
    %T, thickness, resistivity (entered/measured), type, optical constant, calibration,
    %1/64 or 1/1
    for file = 1:length(fileList)
        this_file = fileList{file};
        thick{file,1} = xlsread(this_file,'User','B6');
        res{file,1} = xlsread(this_file,'User','C6');
        oc{file,1} = xlsread(this_file,'User','E6');
        temp{file,1} = 25;
        meas_res{file,1} = xlsread(this_file,'Summary','M2');%'Q2');
        calib{file,1} = xlsread(this_file,'Summary','S2');%'T2');
        doping{file,1} = xlsread(this_file,'Summary','E2');
    end
    info = struct('filename',fileListShort,'thickness',thick,'resistivity',res,'measured_resistivity',meas_res,'optical_constant',oc,'calibration_factor',calib,'temperature',temp,'doping',doping);
    save([dirname '\' samples{index} '\meas_info.mat'],'info');
end

%% Now analyze the data
clear all; close all; clc;
%Process data after HF passivation

dirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 17 2017'; 
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','C-1','C-2','H-1','H-2','FZ'};
samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12'};
lifetime_store = zeros(length(samples),1); 

for i = 1:length(samples)
    load([dirname '\' samples{i} '\Raw_data.mat']); 
    load([dirname '\' samples{i} '\meas_info.mat']); 
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    for j = 1:length(dataSave)
        datanow = dataSave{j}; 
        curves(j)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2); 
        hold all; 
        label{j} = info(j).filename;
        xlim([5e13 1e17])
    end
    xlabel('excess carrier density [cm^-^3]','FontSize',30); 
    ylabel('lifetime [s]','FontSize',30);
    legend(curves',label');
    set(0,'defaultAxesFontSize', 20)
    hgsave(h,[dirname '\' samples{i} '\lifetime summary']);
    print(h,'-dpng','-r0',[dirname '\' samples{i} '\lifetime summary.png']);
    
    %Let's choose measurement 2 and get the lifetime value at 1e15 just for
    %reference
    if length(dataSave)>1
        index = 2; 
    else
        index = 1;
    end
    datanow = dataSave{index}; 
    [deltan,tau] = remove_duplicates(datanow(:,1),datanow(:,2));
    %There's a little bug in this program, for now just do this... 
    try
        lifetime_store(i) = interp1(deltan,tau,1e15); 
    catch
        [deltan,tau] = remove_duplicates(deltan,tau);
        lifetime_store(i) = interp1(deltan,tau,1e15);
    end
end

%% Analyze different states together
clear all; close all; clc;
savedirname = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation';
dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 9 2017';
dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 13 2017'; 
dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 17 2017'; 
dirnames = {dirname1 dirname2 dirname3}; 
% samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','C-1','C-2','H-1','H-2','FZ'};
samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12'};
savename = '_1000s_lifetime summary';
for i = 1:length(samples)
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    curves = [];
    label = {};
    count = 1;
    for k = 1:length(dirnames)
        try 
            load([dirnames{k} '\' samples{i} '\Raw_data.mat']); 
            load([dirnames{k} '\' samples{i} '\meas_info.mat']); 
            flag = 1; 
        catch
            %If that's not successful, just tell me that the sample doesn't
            %exist for that dirname
            warning(['Error accessing data for directory ' num2str(k) ', sample ' samples{i}]);
            flag = 0; 
        end
        if flag == 1
            for j = 1:length(dataSave)
                datanow = dataSave{j}; 
                curves(count)=loglog(datanow(:,1),datanow(:,2),'LineWidth',2); 
                hold all; 
                label{count} = ['Set ' num2str(k) ', #' num2str(j)];
                xlim([5e13 1e17])
                count = count+1; 
            end
        end
    end
    xlabel('excess carrier density [cm^-^3]','FontSize',30); 
    ylabel('lifetime [s]','FontSize',30);
    title(samples{i},'FontSize',30);
    legend(curves',label');
    if strcmp(samples{i},'C-1')==1 || strcmp(samples{i},'C-2')==1 || strcmp(samples{i},'H-1')==1 ||...
        strcmp(samples{i},'H-2')==1 || strcmp(samples{i},'FZ')==1
        ylim([1e-4 1e-2]);
    end
    set(0,'defaultAxesFontSize', 20)
    hgsave(h,[savedirname '\' samples{i} savename]);
    print(h,'-dpng','-r0',[savedirname '\' samples{i} savename '.png']);
end
