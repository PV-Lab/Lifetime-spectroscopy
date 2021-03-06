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
sample_nos = {'19-1','19-2','19-3','19-4','19-5','19-6'...
    '20-1','20-2','20-3','20-4','20-5','20-6',...
    '123-1 1-64','123-2','123-3','123-4 1-1','123-5','123-6',...
    '124-1','124-2','124-3 1-64','124-4','124-5','124-6'};

excel = {'19-1_1-64','19-2_1-64','19-3_1-64','19-4_1-64','19-5_1-64','19-6_1-64'...
    '20-1_1-64','20-2_1-64','20-3_1-64','20-4_1-64','20-5_1-64','20-6_1-64',...
    '123-1_1-64','123-2_1-1','123-3_1-1','123-4_1-1','123-5_1-1','123-6_1-1',...
    '124-1_1-1','124-2_1-1','124-3_1-64','124-4_1-1','124-5_1-1','124-6_1-1'};

lifetime = zeros(length(sample_nos),5); 
filename_start = 'C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Lifetime stage 1\No box\';
for i = 1:length(sample_nos)
    for j = 1:5
        filename = [filename_start sample_nos{i} '\' excel{i} '_' num2str(j) '.xlsm'];
        lifetime(i,j) = xlsread(filename,'Summary','J2'); 
    end
end
