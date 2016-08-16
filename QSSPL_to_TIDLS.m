%This function reads the QSSPL data for a certain sample, where lifetime
%measurements have been taken and processed separately. 
dirname = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160811\A23-6';
%Filenames in order of temperature
filenames = {[dirname '\A23-6_25C.txt'];...
    [dirname '\A23-6_100C.txt'];...
    [dirname '\A23-6_150C.txt'];...
    [dirname '\A23-6_200C.txt']};

%First process into the expected input for TIDLS measurements
%The inputs we need are thickness, doping, type, optical constant,
%temperature, FS used. 

%!!! This is a placeholder
sample_param_filename = 'C:\Users\Mallory\Dropbox (MIT)\TIDLS at UNSW\Advanced system measurements\20160804\For processing\python_inputs.mat';
info = load(sample_param_filename); 

%For each temperature, load the data
for i = 1:length(filenames)
    format_for_TIDLS(filenames{i},sample_param_filename,dirname);
    [tau_mean,deltanq]=average_tau([dirname '\Raw_data.mat'],dirname);
end