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

%Read and plot data from a QSSPL/PC measurement. The files should be .txt
%files with the header removed. 
function [deltanPL,deltanPC,tauPL,tauPC] = read_QSSPL(filename)
    data = textread(filename);
    time = data(:,1);
    deltanPL = data(:,2);
    deltanPC = data(:,2);
    deltanPL = data(:,3);
    tauPC = data(:,4);
    tauPL = data(:,5);
    figure;
    loglog(deltanPL,tauPL,'.');
    hold all;
    loglog(deltanPC,tauPC,'.');
    legend('PL','PC')
    xlabel('Excess carrier density [cm^-^3]','FontSize',20);
    ylabel('Lifetime [s]','FontSize',20); 
end