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