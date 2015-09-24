function plot_lifetime(filename,sample_log,process_sheet,match_sheet)

%This script is written to open Sinton data which has been previously
%processed into a .mat file and parse through the different sets to produce
%plots. 

%load the data of interest
load(filename); 

%How many samples are involved 
[m,n] = size(dataSave); 

%Get the sample log so I know how the samples were processed
[num,txt,raw] = xlsread(sample_log,process_sheet); 
[p,q] = size(raw); 
if q == m+1
    %we probably have a header
    raw(1,:) = []; 
end

%Get the sample log so I know how to match up the samples
[num,txt,raw2] = xlsread(sample_log,match_sheet); 
[p,q] = size(raw2); 
if q == (m/2)+1
    %we probably have a header
    raw2(1,:) = []; 
end

i = 1; 
count = 0; 
indices_complete = []; 
while count<m
    
    data = dataSave{i}; 
    deltan = data(:,1); 
    tau = data(:,2);

    figure('units','normalized','outerposition',[0 0 1 1]) 
    h(1)=loglog(deltan,tau.*1e6,'LineWidth',4); 
    hold all;
    count = count+1; 
    
    %Now we need to find the matching sister
    index = []
    [p,q] = size(raw2);
    file_now = fileListShort{i}; 
    for j = 1:p
        for k = 1:q
            if isempty(strfind(raw2{j,k},file_now(1:5)))==0
                index = [j k];
            end
        end
    end
    
    if isempty(index)==1
        for j = 1:p
            for k = 1:q
                if isempty(strfind(raw2{j,k},file_now(1:4)))==0
                    index = [j k];
                end
            end
        end
    end

    if index(2) ==1
        other_sister = raw2{index(1),2};
    elseif index(2) ==2
        other_sister = raw2{index(1),1}; 
    end
    
    [r,s] = size(raw); 
    index = [];
    for j = 1:r
        if isempty(strfind(raw{j,1},file_now(1:5)))==0
            index = [j];
        end
    end
    
    if isempty(index)==1
        for j = 1:p
            for k = 1:q
                if isempty(strfind(raw{j,1},file_now(1:4)))==0
                    index = [j k];
                end
            end
        end
    end
    
    process1 = raw{index(1),2};
    
    index_sis = [];
    for j = 1:length(fileListShort)
        if isempty(strfind(fileListShort{j},other_sister))==0
            index_sis = j; 
        end
    end
    
    data = dataSave{index_sis}; 
    deltan = data(:,1); 
    tau = data(:,2);
    
    h(2)=loglog(deltan,tau.*1e6,'LineWidth',4); 
    xlabel('Excess carrier density (cm^-^3)','FontSize',30);
    ylabel('Lifetime (\mus)','FontSize',30);
    xlim([5e13 1e17]);
    title([fileListShort{i} ' and ' fileListShort{index_sis}],'FontSize',30); 
    
    index = [];
    for j = 1:r
        if isempty(strfind(raw{j,1},other_sister))==0
            index = j; 
        end
    end
    
    process2 = raw{index,2};
    
    legend(h,process1, process2); 
    count = count+1; 
    
    set(gca,'FontSize',20); 
    set(gca,'LineWidth',2);
    print('-dpng','-r0',[fileListShort{i} '_' fileListShort{index_sis} '.png']);
    
    indices_complete = [indices_complete i index_sis];
    i = i+1; 
    while isempty(find((indices_complete-i)==0))==0
        i =i+1; 
    end
    
    
end

    