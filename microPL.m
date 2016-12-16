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

%Process microPL measurements
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\Larger piece\As brought'; 
[fileList,fileListShort] = getAllFiles(dirname); 
dataStore = cell(size(fileListShort)); 
figure; 
for i = 1:length(fileListShort)
    data = textread(fileList{i}); 
    %Normalize to the starting value
    start_value = max(data(1:10,2)); 
    data(:,2) = data(:,2)./start_value; 
    dataStore{i} = data; 
    plot(data(:,1),data(:,2),'LineWidth',3); 
    hold all; 
end
xlabel('Time','FontSize',20); 
ylabel('Normalized PL counts','FontSize',20); 
legend(fileListShort); 

%Things that I should compare
%Recombination active to non-recombination active grain boundaries
%Intra grain to intra grain regions
%% Compare before/after in situ dark anneal, same location
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation';
before = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\Larger piece\After ex situ dark anneal\broken-7-beforeISDA-532nm-36mW-0P1s-intra-grain-p1.txt';
after = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\Larger piece\After ex situ dark anneal\After in situ dark anneal\broken-7-afterISDA-532nm-36mW-0P1s-intra-grain-p1.txt';
dataStore = cell(2,1); 
data_before = textread(before); 
%Try normalizing to start value
start_value = max(data_before(1:10,2)); 
data_before_norm = data_before;
data_before_norm(:,2) = data_before_norm(:,2)./start_value; 
data_after = textread(after); 
start_value = max(data_after(1:10,2)); 
data_after_norm = data_after;
data_after_norm(:,2) = data_after_norm(:,2)./start_value; 
h=figure;
plot(data_before(:,1),data_before(:,2),'LineWidth',3); 
hold all;
plot(data_after(:,1),data_after(:,2),'LineWidth',3); 
xlabel('Time','FontSize',20); 
ylabel('PL counts','FontSize',20); 
legend('Before','After');
title('Same location before/after in-situ dark anneal','FontSize',20); 
%Plot the normalized data as well
h2=figure;
plot(data_before_norm(:,1),data_before_norm(:,2),'LineWidth',3); 
hold all;
plot(data_after_norm(:,1),data_after_norm(:,2),'LineWidth',3); 
xlabel('Time','FontSize',20); 
ylabel('Normalized PL counts','FontSize',20); 
legend('Before','After');
title('Same location before/after in-situ dark anneal','FontSize',20);
hgsave(h,[dirname '\Insitu_darkanneal']);
print(h,'-dpng','-r0',[dirname '\Insitu_darkanneal.png']);
hgsave(h2,[dirname '\Insitu_darkanneal_norm']);
print(h2,'-dpng','-r0',[dirname '\Insitu_darkanneal_norm.png']);
%% All grain boundaries vs. all intra-grain for in-situ
clear all; close all;
intragrain = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\Intra grain';
grainboundary = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\GB all';
all = figure;
sep_IG = figure;
sep_GB = figure;
all_norm = figure;
sep_IG_norm = figure;
sep_GB_norm = figure; 
%First intragrain 
[fileList,fileListShort] = getAllFiles(intragrain);
dataStore_IG = cell(length(fileListShort),2); 
for i = 1:length(fileListShort)
    data = textread(fileList{i}); 
    dataStore_IG{i,1} = data; 
    figure(all); 
    plot(data(:,1),data(:,2),'k-','LineWidth',2); 
    hold all; 
    figure(sep_IG); 
    plot(data(:,1),data(:,2),'LineWidth',2); 
    hold all; 
    %Normalize to the starting value
    start_value = max(data(1:10,2)); 
    data(:,2) = data(:,2)./start_value; 
    dataStore_IG{i,2} = data;
    figure(all_norm); 
    plot(data(:,1),data(:,2),'k-','LineWidth',2); 
    hold all; 
    figure(sep_IG_norm); 
    plot(data(:,1),data(:,2),'LineWidth',2); 
    hold all; 
end
figure(sep_IG); 
xlabel('Time','FontSize',20); 
ylabel('PL counts','FontSize',20); 
legend(fileListShort); 
hgsave(sep_IG,[intragrain '\IG_only']);
print(sep_IG,'-dpng','-r0',[intragrain '\IG_only.png']);
figure(sep_IG_norm); 
xlabel('Time','FontSize',20); 
ylabel('Normalized PL counts','FontSize',20); 
legend(fileListShort); 
hgsave(sep_IG_norm,[intragrain '\IG_only_norm']);
print(sep_IG_norm,'-dpng','-r0',[intragrain '\IG_only_norm.png']);
%Now grain boundaries
[fileList,fileListShort] = getAllFiles(grainboundary);
dataStore_GB = cell(length(fileListShort),2); 
for i = 1:length(fileListShort)
    data = textread(fileList{i}); 
    dataStore_GB{i,1} = data; 
    figure(all); 
    plot(data(:,1),data(:,2),'b-','LineWidth',2); 
    hold all; 
    figure(sep_GB); 
    plot(data(:,1),data(:,2),'LineWidth',2); 
    hold all; 
    %Normalize to the starting value
    start_value = max(data(1:10,2)); 
    data(:,2) = data(:,2)./start_value; 
    dataStore_GB{i,2} = data;
    figure(all_norm); 
    plot(data(:,1),data(:,2),'b-','LineWidth',2); 
    hold all; 
    figure(sep_GB_norm); 
    plot(data(:,1),data(:,2),'LineWidth',2); 
    hold all; 
end
figure(sep_GB); 
xlabel('Time','FontSize',20); 
ylabel('PL counts','FontSize',20); 
legend(fileListShort); 
hgsave(sep_GB,[grainboundary '\GB_only']);
print(sep_GB,'-dpng','-r0',[grainboundary '\GB_only.png']);
figure(sep_GB_norm); 
xlabel('Time','FontSize',20); 
ylabel('Normalized PL counts','FontSize',20); 
legend(fileListShort); 
hgsave(sep_GB_norm,[grainboundary '\GB_only_norm']);
print(sep_GB_norm,'-dpng','-r0',[grainboundary '\GB_only_norm.png']);
figure(all); 
xlabel('Time','FontSize',20); 
ylabel('PL counts','FontSize',20);
title('Grain boundaries and intragrain together','FontSize',20); 
hgsave(all,[grainboundary '\all']);
print(all,'-dpng','-r0',[grainboundary '\all.png']);
figure(all_norm); 
xlabel('Time','FontSize',20); 
ylabel('Normalized PL counts','FontSize',20);
title('Grain boundaries and intragrain together','FontSize',20); 
hgsave(all_norm,[grainboundary '\all_norm']);
print(all_norm,'-dpng','-r0',[grainboundary '\all_norm.png']);
%% Separating recombination active and non-recombination active grain boundaries (in-situ)
clear all; close all;
recomb = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\GB recombination active'
non_recomb = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\in situ degradation\GB non recombination active';

%% Let's look at the ex-situ measurements
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\ex situ degradation'; 
samples = {'68-3'}; 
%times 
times = {'t0','t1-100s','t2-1000s','t3-10000s'}; 
conditions = {'810nm-3mW-0P5s','810nm-30mW-0P1s'};
points = [4,5;3,3;4,3;2,2];
PL = cell(length(times),length(conditions)); 
for i = 1:length(samples)
    for j = 1:length(times)
        for k = 1:length(conditions)
            PL_store = cell(points(j,k),1);
            x_store = cell(points(j,k),1); 
            h=figure; 
            label = cell(points(j,k),1); 
            for m = 1:points(j,k)
                filename = [dirname '\' samples{i} '-' times{j}...
                    '-' conditions{k} '-p' num2str(m) '.txt']; 
                data = textread(filename); 
                %Throw away the first five points because that's when the
                %scan is starting up
                data(1:5,:) = []; 
                %Center the scan based on the presumed location of the GB
                center = find(data(:,2)==min(data(:,2))); 
                x_center = data(:,1); 
                x_center = x_center(center); 
                %Now adjust the scan based on this center - we want it to
                %be zero
                data(:,1) = data(:,1)-x_center; 
                x_store{m,1} = data(:,1); 
                PL_store{m,1} = data(:,2);
                plot(data(:,1),data(:,2)); 
                label{m} = num2str(m); 
                hold all; 
            end
            x{j,k} = x_store; 
            PL{j,k} = PL_store; 
            xlabel('Distance'); 
            ylabel('PL intensity'); 
            legend(label); 
            title(['Time = ' times{j} ', Condition = ' conditions{k}]);
            %Save the figure
            hgsave(h,[dirname '\' times{j} '_' conditions{k}]);
            print(h,'-dpng','-r0',[dirname '\' times{j} '_' conditions{k} '.png']); 
        end
    end
end
save([dirname '\68-3_exsitu_data.mat'],'x','PL','times','conditions','samples'); 
%% Load the ex-situ data and try averaging
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\ex situ degradation'; 
load([dirname '\68-3_exsitu_data.mat']); 
%We want to plot the different times together
for i = 1:length(conditions)
    h = figure; 
    label = cell(length(times),1); 
    for j = 1:length(times)
        PL_now = PL{j,i}; 
        x_now = x{j,i}; 
        %Now we want to average across each of those points
        [points,nothing] = size(PL_now); 
        x_max = 0; 
        x_min = 0; 
        for k = 1:points
            if x_max < max(x_now{k})
                x_max = max(x_now{k}); 
            end
            if x_min > min(x_now{k})
                x_min = min(x_now{k}); 
            end
        end
        %Now interpolate
        x_interp = linspace(x_min,x_max,500); 
        PL_interp = zeros(500,points); 
        for k = 1:points
            PL_interp(:,k) = interp1(x_now{k},PL_now{k},x_interp); 
        end
        PL_average = nanmean(PL_interp,2); 
        %Now we plot the average
        plot(x_interp,PL_average,'LineWidth',3); 
        hold all; 
        label{j,1} = times{j}; 
    end
    xlabel('Distance'); 
    ylabel('Average PL intensity');  
    title(['Condition = ' conditions{i}]);
    legend(label); 
    %Save the figure
    hgsave(h,[dirname '\Averaged_' conditions{i}]);
    print(h,'-dpng','-r0',[dirname '\Averaged_' conditions{i} '.png']); 
end
%% Look at ex-situ degradation at grains vs. grain boundaries over times
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\ex situ degradation'; 
load([dirname '\68-3_exsitu_data.mat']); 
%We'll want to get the plateau values on either side of the GB + at GB
plateau_rightG = cell(size(PL));
plateau_GB = cell(size(PL)); 
plateau_leftG = cell(size(PL)); 
for i = 1:length(conditions)
    h=figure;
    num_times = [1;100;1000;10000]; 
    for j = 1:length(times)
        PL_now = PL{j,i}; 
        x_now = x{j,i}; 
        [points,nothing] = size(PL_now);
        plateau_leftG_store = zeros(points,1); 
        plateau_GB_store = zeros(points,1); 
        plateau_rightG_store = zeros(points,1); 
        for k = 1:points
            %First pick the area to average over for left GB
            figure;
            plot(x_now{k},PL_now{k}); 
            disp('Select area to average on left-hand side of GB'); 
            [x_range,nothing] = ginput(2); 
            %We take that average and add it to our plateau values
            min_index = find(min(abs(x_now{k}-x_range(1)))==abs(x_now{k}-x_range(1))); 
            max_index = find(min(abs(x_now{k}-x_range(2)))==abs(x_now{k}-x_range(2)));
            PL_plateau = PL_now{k}; 
            plateau_leftG_store(k,1) = nanmean(PL_plateau(min_index:max_index,1));
            %Now we repeat for the right plateau
            disp('Select area to average on right-hand side of GB'); 
            [x_range,nothing] = ginput(2);
            %We take that average and add it to our plateau values
            min_index = find(min(abs(x_now{k}-x_range(1)))==abs(x_now{k}-x_range(1))); 
            max_index = find(min(abs(x_now{k}-x_range(2)))==abs(x_now{k}-x_range(2)));
            plateau_rightG_store(k,1) = nanmean(PL_plateau(min_index:max_index,1));
            %For the grain boundary we just find the minimum value
            plateau_GB_store(k,1) = min(PL_plateau); 
        end
        figure(h); 
        plateau_rightG{j,i} = plateau_rightG_store; 
        plateau_GB{j,i} = plateau_GB_store; 
        plateau_leftG{j,i} = plateau_leftG_store; 
        %We also want to plot these together
        leg(1)=semilogx(num_times(j).*ones(size(plateau_leftG_store)),plateau_leftG_store,'ro','MarkerSize',10,'LineWidth',3); 
        hold all; 
        leg(2)=semilogx(num_times(j).*ones(size(plateau_rightG_store)),plateau_rightG_store,'go','MarkerSize',10,'LineWidth',3); 
        hold all; 
        leg(3)=semilogx(num_times(j).*ones(size(plateau_GB_store)),plateau_GB_store,'bx','MarkerSize',10,'LineWidth',3);  
        hold all; 
    end
    xlabel('Times [s]','FontSize',20); 
    ylabel('PL intensity [counts]','FontSize',20); 
    legend(leg,{'Left grain';'Right grain';'Grain boundary'}); 
    title(['Condition = ' conditions{i}]);
    %Save the figure
    hgsave(h,[dirname '\PlateauValues_' conditions{i}]);
    print(h,'-dpng','-r0',[dirname '\PlateauValues_' conditions{i} '.png']);
end
save([dirname '\PlateauValues.mat'],'plateau_rightG','plateau_GB','plateau_leftG','num_times','times','conditions');         
%% Just plot the plateau values from before
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\Mallory in Australia\PERC LeTID\microPL!\text data\ex situ degradation'; 
load([dirname '\PlateauValues.mat']); 
for i = 1:length(conditions)
    h=figure;
    for j = 1:length(times)
        plateau_rightG_store = plateau_rightG{j,i};
        plateau_GB_store = plateau_GB{j,i}; 
        plateau_leftG_store = plateau_leftG{j,i};
        %We also want to plot these together
        leg(1)=semilogx(num_times(j).*ones(size(plateau_leftG_store)),plateau_leftG_store,'ro','MarkerSize',10,'LineWidth',3); 
        hold all; 
        leg(2)=semilogx(num_times(j).*ones(size(plateau_rightG_store)),plateau_rightG_store,'go','MarkerSize',10,'LineWidth',3); 
        hold all; 
        leg(3)=semilogx(num_times(j).*ones(size(plateau_GB_store)),plateau_GB_store,'bx','MarkerSize',10,'LineWidth',3);  
        hold all; 
    end
    xlabel('Times [s]','FontSize',20); 
    ylabel('PL intensity [counts]','FontSize',20); 
    legend(leg,{'Left grain';'Right grain';'Grain boundary'}); 
    title(['Condition = ' conditions{i}]);
    %Save the figure
    hgsave(h,[dirname '\PlateauValues_' conditions{i}]);
    print(h,'-dpng','-r0',[dirname '\PlateauValues_' conditions{i} '.png']);
end