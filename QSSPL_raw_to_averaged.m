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

%% Average raw data from QSSPL system - put in a format that we can read in the QSSPL analyzer 
clear all; close all; clc;

%We need to enter the inputs such that these match our filenames
%Check the filename structure in average_QSSPL first

%This is the directory where our measurements are stored
dirname = 'C:\Users\Mallory Jensen\Documents\Noncontact crucible\TIDLS UNSW\Data\sorted by sample\16 6 28 N 1\25C\9-26-16';

%This is the sample name in the filename, dirname\sample name
sample = '16-6-28-N-1'; 

%This is the part that comes between the temperature value (if applicable)
%and the rest of the filename
Ctonum = 'glass_1-1_';

%This is the time measured before the flash arrives, used to "zero" the
%data
time_before = .001; %s

%These are the temperatures at which the measurements are taken. All
%measurements must be in the same directory for this to work. 
T = [25]; 

%This is the number of measurements taken under that condition. The size of
%this vector should match T - loop over both of these together. 
num_repeats = [5]; 

%With the above inputs, call the function that zeros and averages the data
[dataSave] = average_QSSPL(dirname,sample,T,num_repeats,time_before,Ctonum);

%Now, given this data, we want to write this to a new file which can be
%read into QSSPL
%Note that you may get an error here if the filename is unusual. You want
%the first input of copyfile in write_QSSPL_analyser (lines 36-37) to match
%what the structure was in average_QSSPL
for i = 1:length(T)
    matrix_to_write = dataSave{i}; 
    write_QSSPL_analyser(dirname,sample,T(i),matrix_to_write,Ctonum);
end
    
