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

%%Revise the thickness, resistivity, and optical constant in Sinton
%%spreadsheets
clear all; close all; clc; 
meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Repassivated samples\Degradation\measurement_details.xlsx'; 
%Get the measurement details
[times,filenames] = xlsread(meas_details,'filenames'); 
[meas,samples] = xlsread(meas_details,'measurements');
samples(1,:) = [];
%Make these the same size
filenames = filenames(2:end,2);

%Get the sample parameters
sample_params = 'C:\Users\Mallory Jensen\Documents\LeTID\Dartboard\Sample measurements.xlsx'; 
[params,names_params] = xlsread(sample_params,'sample summary','A2:D13'); 

for i = 1:length(samples)
    if strcmp(samples{i},'C-h-5')==0 && strcmp(samples{i},'C-L-5')==0
        meas_thissample = meas(i,:);
        %Get the proper parameters
        params_index = find(strcmp(names_params,samples{i})==1); 
        params_this_sample = params(params_index,:); 
        for j = 1:length(meas_thissample)
            if isnan(meas_thissample(j))==0
                %now create the proper directory name
                findex = find(meas_thissample(j)==times);  
                dirname = [filenames{findex} '\' samples{i}];
                %Get all of the files in this directory
                fileList = dir(fullfile(dirname,'*.xlsm')); 
                fileList = {fileList.name}'; 
                for k = 1:length(fileList)
                    filename = fullfile(dirname,fileList{k}); 
                    try
                        %Write the thickness
                        xlswrite(filename,params_this_sample(1),'User','B6'); 
                    catch
                        disp(['There was an error with ' filename ' thickness']); 
                    end
                    try
                        %Write the optical constant
                        xlswrite(filename,params_this_sample(3),'User','E6'); 
                    catch
                        disp(['There was an error with ' filename ' optical constant']); 
                    end
                    try
                        %Write the resistivity
                        xlswrite(filename,params_this_sample(2),'User','C6'); 
                    catch
                        disp(['There was an error with ' filename ' resistivity']); 
                    end
                end
            end
        end
    end
end
