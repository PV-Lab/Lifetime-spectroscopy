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

clear all; close all; clc; 
%Run through different temperatures and evaluate optical constant
dirname = 'C:\Users\Mallory\Documents\Australia\Optical constant\';
%folder names
folders = {'25C','75C','100C','175C','200C'}; 
T = [25 75 100 175 200]; 
OC_range = [0.8 1.2]; 
resolution = 20; 
inj_range = [3e15 1e16]; %assume this can be the same for all temperatures
for i = 1:length(folders)
    location = [dirname folders{i}];
    [fileList,fileListShort] = getAllFiles(location);
    %gen is first (1-1)
    filename_gen = fileList{1}; 
    filename_trans = fileList{2}; 
    [deltan_trans,tau_trans,dataSave_gen,OC_values,MSE,min_index]=modify_OC(filename_trans,filename_gen,OC_range,inj_range,resolution);
    OC_win(i) = OC_values(min_index); 
end

figure; 
plot(T,OC_win,'o'); 
xlabel('Temperature [C]','FontSize',20); 
ylabel('Best optical constant [-]','FontSize',20); 
title('Float zone test sample','FontSize',20); 
