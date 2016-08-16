%This function reads the QSSPL data (raw) for a certain sample, where lifetime
%measurements have been taken and processed separately.
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160811\A23-6';
%Filenames in order of temperature
filenames = {[dirname '\A23-6_25C.txt'];...
    [dirname '\A23-6_100C.txt'];...
    [dirname '\A23-6_150C.txt'];...
    [dirname '\A23-6_200C.txt']};
%Sample parameters
% temperature = {};
% thickness_samp = ;
% doping_samp = ;
% type_samp = 'p'; 
% OC = {0.7;0.7;0.7;0.7};
% FS_samp = ;
% 
% %First process into the expected input for TIDLS measurements
% %The inputs we need are thickness, doping, type, optical constant,
% %temperature, FS used.
% %type, FS, thickness, doping are all the same regardless of T
% type = cell(size(filenames)); FS = cell(size(filenames)); 
% thickness = cell(size(filenames)); doping = cell(size(filenames));
% for i = 1:length(filenames)
%     type{i} = type_samp;
%     FS{i} = FS_samp;
%     thickness{i} = thickness_samp;
%     doping{i} = doping_samp;
% end

%!!! This is a placeholder
sample_param_filename = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160804\For processing\python_inputs.mat';
info = load(sample_param_filename); 

%For each temperature, load the data
for i = 1:length(filenames)
    format_for_TIDLS(filenames{i},sample_param_filename,dirname);
    [tau_mean,deltanq]=average_tau([dirname '\Raw_data.mat'],dirname);
end

%% As an alternative, try averaging the raw data. But first we need to figure out how to correct the raw data
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160811\A23-6';
file = [dirname '\A23-6_25C_1.Raw Data.dat'];
savename = [dirname '\test.mat'];
[data_now,time,PC,PL,gen]=format_QSSPL_rawdata(file,savename);

%Define the wait time before the pulse arrives
time_before = .001; %s
[PC_corr] = correct_dark_voltage(time,PC,time_before);
% %Let's take 75% of that just to be conservative
% time_before = time_before*0.75; 
% %Now look for the points to find the dark voltage on PC
% PC_dark = PC(find(time<=time_before)); 
% PC_dark = nanmean(PC_dark); 
% %Now correct all values for the dark voltage
% PC = PC-PC_dark; 
% %Now replot the results
% figure;
% plot(time,gen);
% hold all;
% plot(time,PC);
% hold all;
% plot(time,PL);
% xlabel('Time [s]','FontSize',20); 
% ylabel('Voltage [V]','FontSize',20);
% legend('Generation','PC','PL'); 
% title('Dark voltage corrected','FontSize',20); 
%% Now try averaging for a series 
clear all; close all; 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160811\A23-6';
sample = 'A23-6'; 
time_before = .001; %s
T = [25 100 150 200]; 
num_repeats = [5 10 10 10]; 
[dataSave] = average_QSSPL(dirname,sample,T,num_repeats,time_before);

%Now, given this data, we want to write this to a new file which can be
%read into QSSPL


% for i = 1:length(T)
%     for j =1:num_repeats(i)
%         %Make filename
%         filename = [dirname '\' sample '_' num2str(T(i)) 'C_' num2str(j)...
%             '.Raw Data.dat'];
%         %Make the savename
%         savename = [dirname '\' sample '_' num2str(T(i)) 'C_' num2str(j)...
%             '_rawData.mat'];
%         %Read the data
%         [data_now,time,PC,PL,gen]=format_QSSPL_rawdata(filename,savename);
%         %Correct the PC
%         [PC_corr] = correct_dark_voltage(time,PC,time_before);
%         %Let's assume the time vector is the same for all of the
%         %measurements
%         if j == 1
%             PC_sum = PC_corr; 
%             gen_sum = gen; 
%             time_store = time; 
%         else
%             PC_sum = PC_sum+PC_corr; 
%             gen_sum = gen_sum+gen; 
%         end
%         if isequal(time,time_store)==0
%             disp('Error with time vectors, times are not equivalent sample-to-sample.'); 
%         else
%             time_store = time; 
%         end
%     end
%     %Now we have the average for our given temperature
%     PC_overall = PC_sum./num_repeats(i);
%     gen_overall = gen_sum./num_repeats(i); 
%     %Plot the result
%     figure;
%     plot(time_store,gen_overall);
%     hold all;
%     plot(time_store,PC_overall);
%     xlabel('Time [s]','FontSize',20); 
%     ylabel('Voltage [V]','FontSize',20);
%     legend('Generation','PC'); 
%     title(['Average PC, T = ' num2str(T(i)) 'C'],'FontSize',20); 
% end
%             
        
