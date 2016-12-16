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

function [possible_MCD,Joe]=Joe_curve(filename,min_MCD,max_MCD,doping)
%This function takes as input the filename (+location) of a new Sinton
%QSSPC spreadsheet and the minimum and maximum carrier densities. The
%function obtains a curve of Joe by manipulating the Excel spreadsheet -
%writing MCD values and then reading the result. The last input is the
%doping level which is used for plotting a reference value for the Joe
%curves. 
%Define the range of possible MCDs
possible_MCD = linspace(min_MCD,max_MCD,40);
for i = 1:length(possible_MCD)
    %Write the MCD to the file
    xlswrite(filename,possible_MCD(i),'User','F6');
    %Now that we've written it, go in a read the Joe
    Joe(i) = xlsread(filename,'Summary','L2'); 
end
%Plot the results
figure;
plot(possible_MCD,Joe,'o');
hold all; 
plot([doping doping],[0 max(Joe)],'k-');
xlabel('Excess carrier density [cm^-^3]','FontSize',20);
ylabel('J_o_e [A/cm^2]','FontSize',20); 
