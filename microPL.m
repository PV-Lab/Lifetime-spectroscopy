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
%Before/after in situ dark anneal, same location