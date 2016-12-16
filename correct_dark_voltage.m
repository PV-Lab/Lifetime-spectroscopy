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

function [PC_corr] = correct_dark_voltage(time,PC,time_before)
time_before = time_before*0.75; 
%Now look for the points to find the dark voltage on PC
PC_dark = PC(find(time<=time_before)); 
PC_dark = nanmean(PC_dark); 
%Now correct all values for the dark voltage
PC_corr = PC-PC_dark; 
%Now replot the results
figure;
plot(time,PC);
hold all;
plot(time,PC_corr);
xlabel('Time [s]','FontSize',20); 
ylabel('Voltage [V]','FontSize',20);
legend('PC','PC corrected'); 
title('Dark voltage corrected','FontSize',20);