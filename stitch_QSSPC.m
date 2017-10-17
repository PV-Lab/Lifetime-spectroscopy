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

%This function takes as input a directory with QSSPC files in it (same
%sample, different measurement "type") and then stitches the data from
%those files together into one .mat file

function [dataSave] = stitch_QSSPC(dirname)

%Get all of the files in that directory with extension .xlsm
fileList = dir(fullfile(dirname,'*.xlsm')); 
fileList = {fileList.name}';

process_xls_data(dirname,[dirname '\Raw_data_before_stitch.mat']);

%Load the data we just saved
load([dirname '\Raw_data_before_stitch.mat']); 
deltan = []; 
tau = []; 

%We also need the relevant data. Check as we go to see if this is the same for
%each measurement
for file = 1:length(fileList)
    %Put the relevant data in the vectors
    data_now = dataSave{file}; 
    deltan = [deltan; data_now(:,1)]; 
    tau = [tau; data_now(:,2)]; 
    this_file = fullfile([dirname],fileList{file});
    thick{file,1} = xlsread(this_file,'User','B6');
    res{file,1} = xlsread(this_file,'User','C6');
    oc{file,1} = xlsread(this_file,'User','E6');
    temp{file,1} = 25;
    meas_res{file,1} = xlsread(this_file,'Summary','N2');
    calib{file,1} = xlsread(this_file,'Summary','T2');
    doping{file,1} = xlsread(this_file,'Summary','E2');
    if file > 1
        if thick{file,1}~=thick{1,1} || res{file,1}~=res{1,1} || ...
                oc{file,1}~=oc{1,1} || doping{file,1}~= doping{1,1}
            %Tell the user we hav an issue
            disp('Possible error: The sample parameters are not consistent within this directory'); 
        end
    end
end

%Replace dataSave 
dataSave = [deltan,tau]; 
%Now sort so we don't plot weird stuff
dataSave = sort(dataSave,1,'ascend'); 
%Plot the result
h=figure('units','normalized','outerposition',[0 0 1 1]);
loglog(dataSave(:,1),dataSave(:,2),'LineWidth',2); 
xlabel('Excess carrier density (cm^-^3)','FontSize',20);
ylabel('Lifetime (seconds)','FontSize',20);
%DataSave should be a cell so we can use it like we would in other scripts
dataSave = {dataSave}; 

%Now save everything
info = struct('filename',fileList{1},'thickness',thick{1},'resistivity',...
    res{1},'measured_resistivity',meas_res{1},'optical_constant',oc{1},...
    'calibration_factor',calib{1},'temperature',temp{1},'doping',doping{1});
save([dirname '\meas_info.mat'],'info');
fileListShort = [dirname '\' fileList{1}]; 
save([dirname '\Raw_data.mat'],'dataSave','fileListShort'); 
end