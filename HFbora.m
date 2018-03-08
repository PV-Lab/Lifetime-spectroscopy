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
    plotting_groups,plotting_names,meas_details,max_time,...
    latesta,latestb]=HFbora(bora)
%This function takes as input a string, 'set-a' or 'set-b' or 'compE' defining the
%sample set to be analyzed. This is specific to a particular experiment.
%Then we define the directory names that we want to plot together on a net
%lifetime curve. This should be updated whenever we have new data we want
%to analyze. 

%First, sample numbers, then dirnames to plot
if strcmp(bora,'set-a')==1
    meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_seta.xlsx'; 
    samples = {'44a','45a','49a','50a','52a','53a','54a','55a','56a','60a','61a','H-1','H-2','FZ','FZ-12','66-2','FZ-new','FZ-new2'}; 
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
    dirname43 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 24 2017';
    dirname44 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 26 2017';
    dirname45 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 31 2017';
    dirname46 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 2 2017';
    dirname47 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 7 2017';
    dirname48 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 9 2017';
    dirname49 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 8 2018';
    dirname50 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 11 2018';
    dirname51 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 19 2018';
    dirname52 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 23 2018';
    dirnames = {dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 dirname8 ...
        dirname10 dirname11 dirname12 dirname13 dirname14 dirname15 dirname16 ...
        dirname17 dirname18 dirname19 dirname20 dirname21 dirname22 dirname23 ...
        dirname24 dirname25 dirname26 dirname27 dirname28 dirname29 dirname30 ...
        dirname31 dirname32 dirname33 dirname34,dirname35,dirname36 dirname37 ...
        dirname38 dirname39 dirname40 dirname41  dirname42 dirname43 ...
        dirname44 dirname45 dirname46 dirname47 dirname48 dirname49 dirname50 ...
        dirname51 dirname52}; 
    labels = {'initial','1000s','2000s','3000s','4000s','5000s','10000s',...
        '20000s','30000s','40000s','50000s','60000s','70000s','80000s','90000s',...
        '100000s','154495','206005','258325','300025','349990','402730s',...
        '454750s','508360s' '601540s', '702100s','801610s','902260s',...
        '1004410s','1167760s','1408600s','1579630s','1804090s','2024170s',...
        '2509510s','2536210s','2561230s','3140860s' '3140861s','3196270s',...
        '3252190s','3308020s','3370630s','3428710s','3488410s','3540430s',...
        '3540430s check','3612910s','3687370s','3761830s'};
    savename = '_seta_3761830s';
    surface_control = {'FZ','FZ-new','FZ-new2'};
    control = {'H-1','H-2','FZ','FZ-12','68-2','66-2','FZ-new','FZ-new2';...
    'Unfired Cz (120 min H)','Fired Cz (120 min H)','FZ passivation',...
    'FZ degradation','mc-Si 950C fired undegraded',...
    'mc-Si 750C fired undegraded','new FZ passivation','FZ passivation #3'};
    fired = {'49a','53a','56a','52a','55a','60a';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    unfired = {'61a','54a','50a','45a','44a','60a';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    plotting_groups = {control, fired, unfired}; 
    plotting_names = {'controls','fired','unfired'}; 
    max_time = 3761830; 
    latesta = []; latestb = [];
elseif strcmp(bora,'set-b')==1
    meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_setb.xlsx'; 
    samples = {'44b','45b','49b','50b','52b','53b','54b','55b','56b','60b','61b','C-1','C-2','66-2','FZ-new','FZ','FZ-12','FZ-new2'};
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
    dirname13 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 24 2017'; 
    dirname14 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 26 2017'; 
    dirname15 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 31 2017';
    dirname16 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 2 2017';
    dirname17 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 7 2017';
    dirname18 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 9 2017';
    dirname19 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 8 2018';
    dirname20 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 11 2018';
    dirname21 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 19 2018';
    dirname22 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 23 2018';
    dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
        dirname8 dirname9 dirname10 dirname11 dirname12 dirname13 dirname14 ...
        dirname15 dirname16 dirname17 dirname18 dirname19 dirname20 ...
        dirname21 dirname22}; 
    labels = {'initial' '50s' '100s' '250s' '500s' '501s' '750s' '1000s' '1500s' ...
        '2000s' '2500s' '3000s' '3500s' '4000s' '4500s' '5000s' '6000s' ...
        '7000s' '7000s check' '10000s' '60280s' '113020s'};
    savename = '_setb_113020s';
    surface_control = {'FZ','FZ-new','FZ-new2'}; 
    control = {'C-1','C-2','FZ','66-2','FZ-12','FZ-new','FZ-new2';...
    'Unfired Cz (no H)','Fired Cz (no H)','FZ passivation','mc-Si 750C fired undegraded','FZ degradation','FZ passivation #2','FZ passivation #3'};
    fired = {'50b','54b','44b','49b','61b','56b';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    unfired = {'52b','45b','55b','60b','53b','56b';...
        '0 min','10 min','30 min','120 min','30 min no H','LeTID control'};
    plotting_groups = {control, fired, unfired}; 
    plotting_names = {'controls','fired','unfired'}; 
    max_time = 113020; 
    latesta = []; latestb = [];
elseif strcmp(bora,'compE')==1
    meas_details = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\measurement_details_compE.xlsx'; 
    samples = {'66-2','FZ-new','68-2','68-4','FZ-12','FZ','FZ-new2'};
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
    dirname23 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 24 2017';
    dirname24 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 26 2017';
    dirname25 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\October 31 2017';
    dirname26 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 2 2017';
    dirname27 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 7 2017';
    dirname28 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\November 9 2017';
    dirname29 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 8 2018';
    dirname30 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 11 2018';
    dirname31 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 19 2018';
    dirname32 = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\January 23 2018';
    dirnames = {dirname1 dirname2 dirname3 dirname4 dirname5 dirname6 dirname7 ...
        dirname8 dirname9 dirname10 dirname11 dirname12 dirname13 dirname14 ...
        dirname15 dirname16 dirname17 dirname18 dirname19 dirname20 dirname21 ...
        dirname22 dirname23 dirname24 dirname25 dirname26 dirname27 dirname28 ...
        dirname29 dirname30 dirname31 dirname32};
    labels = {'initial','81330s','103800s','129780s','151350s','187140s',...
        '219330s','247410s','273060s','299760s','324780s','376890s',...
        '434940s','485940s','543570s','599160s','661830s','724170s',...
        '781890s','839070s','897840s','953250s','1009170s','1065000s',...
        '1127610s','1185720s','1245420s','1297440s','1297440s check',...
        '1369920s','1444380s','1518840s'};
    plotting_groups = {{'FZ-new','FZ-12','FZ','FZ-new2';...
        'FZ passivation new','FZ degradation','FZ passivation','FZ passivation #3'},...
        {'66-2','68-2','68-4';'mc-Si 750C undegraded',...
        'mc-Si 950C bare','mc-Si 950C passivated'}}; 
    plotting_names = {'controls','compE'}; 
    surface_control = {'FZ-new','FZ','FZ-new2'};
    savename = '_compE_1518840s';
    max_time = 1518840; 
    latesta = []; latestb = [];
elseif strcmp(bora,'compare')==1
    %Latest files
    latesta = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set a\3761830s\set-a_all_data_seta_3761830s.mat';
    latestb = 'C:\Users\Mallory Jensen\Documents\LeTID\Hydrogenation experiment\HF passivation\set b\113020s\set-b_all_data_setb_113020s.mat'; 
    %This is just for the degradation comparison, not the lifetime plotting
    plotting_groups = {{'61a','52b';'a','b'},{'54a','45b';'a','b'},...
        {'50a','55b';'a','b'},{'45a','60b';'a','b'},{'44a','53b';'a','b'},...
        {'49a','50b';'a','b'},{'53a','54b';'a','b'},{'56a','44b';'a','b'},...
        {'52a','49b';'a','b'},{'55a','61b';'a','b'},{'60a','56b';'a','b'},...
        {'C-1','H-1';'b','a'},{'C-2','H-2';'b','a'}};
    plotting_names = {'unfired 0 min','unfired 10 min','unfired 30 min',...
        'unfired 120 min','unfired 30 min no H','fired 0 min','fired 10 min',...
        'fired 30 min','fired 120 min',' fired 30 min no H','LeTID control',...
        'Cz unfired','Cz fired'};
    savename = '_comparison_maxTimes';
    max_time = 150000; %just show what will be comparable between the two
    %These outputs don't matter for this case
    samples={};dirnames={};labels={};surface_control={};
    meas_details=[];
else
    disp('There was an error assigning the sample set'); 
end
