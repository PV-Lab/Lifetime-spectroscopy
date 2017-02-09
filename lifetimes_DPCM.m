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

%Make lifetime values for DPCM;
 load('C:\Users\Mallory Jensen\Documents\LeTID\Experiment 0\rev resistivity, OC\processed_data_8e14_20161122.mat');
 
 %which sample?
 sample_index = 2; 
 
 %which times?
 time_index = [15 47 74 79]; %intermediate, max deg, initial regen, later regen
 
 lifetime_store = cell(size(time_index)); 
 deltan_store = cell(size(time_index)); 
 
 figure;
 
 for i = 1:length(time_index)
     datanow = dataSave{sample_index,time_index(i)}; 
     h(i)=loglog(datanow(:,1),datanow(:,2)); 
     hold all; 
     %Let's sample the data
     min_inj = min(datanow(:,1)); 
     if min_inj<1e14
         min_inj = 1e14; 
     end
     max_inj = 2e16; 
     %Find where those are
     mindiff = abs(datanow(:,1)-min_inj);
     min_index = find(mindiff==min(mindiff)); 
     maxdiff = abs(datanow(:,1)-max_inj); 
     max_index = find(maxdiff==min(maxdiff)); 
     abs_min = datanow(:,1); abs_min = abs_min(min_index); 
     abs_max = datanow(:,1); abs_max = abs_max(max_index); 
     deltans = logspace(log10(abs_min),log10(abs_max),10); 
     lifetimes = interp1(datanow(:,1),datanow(:,2),deltans);
     if min_inj==abs_min
         %then this is truly the first index and we should have a nan
         if isnan(lifetimes(1))==1
             lifetime_now = datanow(:,2);
             lifetimes(1)=lifetime_now(1); 
         end
     end
     lifetime_store{i} = lifetimes';
     deltan_store{i} = deltans';
     hold all;
     loglog(deltans,lifetimes,'o'); 
 end
 
 xlabel('excess carrier density [cm^-^3]'); 
 ylabel('lifetime [s]'); 
 legend(h,num2str(time_index')); 
 