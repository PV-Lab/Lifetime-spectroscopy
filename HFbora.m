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
%This script performs analysis for samples degraded and measured with HF
%passivation. This includes improvements for switching between different
%sample sets as well as for analyzing the lifetime independent of surface
%fluctuations. 

function [samples,dirnames,labels,savename,surface_control,...
    plotting_groups,plotting_names,meas_details,max_time]=HFbora(bora)
%This function takes as input a string, 'set-a' or 'set-b' or 'compE' defining the
%sample set to be analyzed. This is specific to a particular experiment.
%Then we define the directory names that we want to plot together on a net
%lifetime curve. This should be updated whenever we have new data we want
%to analyze. 

%First, sample numbers, then dirnames to plot
if strcmp(bora,'set-a')==1
    meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_seta.xlsx'; 
    samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12','66-2'}; 
    dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 9 2017';
    dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 13 2017'; 
    dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 17 2017'; 
    dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 19 2017';
    dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 23 2017';
    dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 25 2017';
    dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 27 2017';
    dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 31 2017';
    dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 2 2017';
    dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 6 2017';
    dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 8 2017';
    dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 13 2017';
    dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 15 2017';
    dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\February 17 2017';
    dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 1 2017\revised calibration';
    dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 8 2017';
    dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 11 2017';
    dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 15 2017';
    dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 21 2017';
    dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 23 2017';
    dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\March 30 2017';
    dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 3 2017';
    dirname23 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 11 2017';
    dirname24 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 13 2017';
    dirname25 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 20 2017';
    dirname26 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 21 2017';
    dirname27 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\April 27 2017';
    dirname28 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 3 2017';
    dirname29 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 5 2017';
    dirname30 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 9 2017';
    dirname31 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 12 2017';
    dirname32 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 19 2017';
    dirname33 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 25 2017';
    dirname34 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 1 2017';
    dirname35 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 7 2017';
    dirname36 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 13 2017';
    dirname37 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\July 12 2017';
    dirname38 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017';
    dirname39 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 11 2017';
    dirname40 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 14 2017';
    dirname41 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 18 2017';
    dirname42 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 19 2017';
    dirnames = {dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 dirname8 ...
        dirname10 dirname11 dirname12 dirname13 dirname14 dirname15 dirname16 ...
        dirname17 dirname18 dirname19 dirname20 dirname21 dirname22 dirname23 ...
        dirname24 dirname25 dirname26 dirname27 dirname28 dirname29 dirname30 ...
        dirname31 dirname32 dirname33 dirname34,dirname35,dirname36 dirname37 ...
        dirname38 dirname39 dirname40 dirname41 dirname42}; 
    labels = {'initial','1000s','2000s','3000s','4000s','5000s','10000s',...
        '20000s','30000s','40000s','50000s','60000s','70000s','80000s','90000s',...
        '100000s','154495','206005','258325','300025','349990','402730s',...
        '454750s','508360s' '601540s', '702100s','801610s','902260s',...
        '1004410s','1167760s','1408600s','1579630s','1804090s','2024170s',...
        '2509510s','2536210s','2561230s','3140860s' '3140861s','3196270s'};
    savename = '_seta_3196270s';
    surface_control = {'FZ','FZ-new'};; 
    control = {'H-1','H-2','FZ','FZ-12','68-2','66-2','FZ-new';...
    'Unfired Cz (120 min H)','Fired Cz (120 min H)','FZ passivation',...
    'FZ degradation','mc-Si 950C fired undegraded',...
    'mc-Si 750C fired undegraded','new FZ passivation'};
    fired = {'49a','53a','56a','52a','55a','60a';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    unfired = {'61a','54a','50a','45a','44a','60a';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    plotting_groups = {control, fired, unfired}; 
    plotting_names = {'controls','fired','unfired'}; 
    max_time = 3196270; 
elseif strcmp(bora,'set-b')==1
    meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_setb.xlsx'; 
    samples = {'44b','45b','49b','50b','52b','53b','54b','55b','56b','60b','61b','C-1','C-2','66-2','FZ-new','FZ','FZ-12'};
    dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017'; 
    dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 11 2017';
    dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 16 2017';
    dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 29 2017';
    dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 30 2017';
    dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 8 2017';
    dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 14 2017';
    dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 20 2017';
    dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 22 2017';
    dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 28 2017';
    dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 12 2017';
    dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 19 2017';
    dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
        dirname8 dirname9 dirname10 dirname11 dirname12}; 
    labels = {'initial' '50s' '100s' '250s' '500s' '501s' '750s' '1000s' '1500s' ...
        '2000s' '2500s' '3000s'};
    savename = '_setb_3000s';
    surface_control = {'FZ','FZ-new'}; 
    control = {'C-1','C-2','FZ','66-2','FZ-12','FZ-new';...
    'Unfired Cz (no H)','Fired Cz (no H)','FZ passivation','mc-Si 750C fired undegraded','FZ degradation','FZ passivation #2'};
    fired = {'50b','54b','44b','49b','61b','56b';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    unfired = {'52b','45b','55b','60b','53b','56b';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    plotting_groups = {control, fired, unfired}; 
    plotting_names = {'controls','fired','unfired'}; 
    max_time = 3000; 
elseif strcmp(bora,'compE')==1
    meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_compE.xlsx'; 
    samples = {'66-2','FZ-new','68-2','68-4','FZ-12','FZ'};
    dirname1 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 5 2017';
    dirname2 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 9 2017';
    dirname3 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 12 2017';
    dirname4 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 19 2017';
    dirname5 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\May 25 2017';
    dirname6 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 1 2017';
    dirname7 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 7 2017';
    dirname8 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\June 13 2017';
    dirname9 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\July 12 2017';
    dirname10 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 9 2017';
    dirname11 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 11 2017';
    dirname12 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 14 2017';
    dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 16 2017';
    dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 29 2017';
    dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\August 30 2017';
    dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 8 2017';
    dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 14 2017';
    dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 20 2017';
    dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 22 2017';
    dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\September 28 2017';
    dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 12 2017';
    dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 19 2017';
    dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
        dirname8 dirname9 dirname10 dirname11 dirname12 dirname13 dirname14 ...
        dirname15 dirname16 dirname17 dirname18 dirname19 dirname20 dirname21 ...
        dirname22};
    labels = {'initial','81330s','103800s','129780s','151350s','187140s',...
        '219330s','247410s','273060s','299760s','324780s','376890s',...
        '434940s','485940s','543570s','599160s','661830s','724170s',...
        '781890s','839070s','897840s','953250s'};
    plotting_groups = {samples}; 
    plotting_names = {'compE'}; 
    surface_control = {'FZ-new','FZ'};
    savename = '_compE_953250s';
    max_time = 953250; 
else
    disp('There was an error assinging the sample set'); 
end
