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

function [k] = calculate_Mod(T)
%This function calculates the k-value for a molybdenum donor defect at
%a given temperature (K). The fit is based on equations 6 and 7 in Paudyal JAP
%2010. Input can be a single temperature value or a vector of temperatures.
%The output structure will correspond to the input structure. 

for i = 1:length(T)
    k(i) = (6.9587e5)*(T(i)^-1.88); 
end

figure;
semilogy(T,k,'o','LineWidth',3); 
xlabel('Temperature [K]','FontSize',20); 
ylabel('k-value [-]','FontSize',20);
title('Mo_d','FontSize',20); 