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

%This function writes data to .dat and .inf files which are compatible with
%QSSPL analyzer. 
function write_QSSPL_analyser(dirname,sample,T,matrix_to_write,Ctonum)
%     filename_write = [dirname '\' sample '_' num2str(T) Ctonum 'averaged.Raw Data.dat'];
    filename_write = [dirname '\' sample '_' num2str(T) 'C' Ctonum 'averaged.Raw Data.dat'];
    %We need to write the first line
    first_line = '# Time (s)	Generation (V)	PC (V)	PL (V)';
    dlmwrite(filename_write,first_line,'delimiter',''); 
    %Now append the data
    dlmwrite(filename_write,matrix_to_write,'-append','delimiter','\t','precision','%-6.18e','newline','pc'); 
    %Make the .inf file
%     copyfile([dirname '\' sample '_' num2str(T) Ctonum '1.inf'],[dirname '\' sample '_' num2str(T) Ctonum 'averaged.inf']);
    copyfile([dirname '\' sample '_' Ctonum '1.inf'],[dirname '\' sample '_' num2str(T) 'C' Ctonum 'averaged.inf']);
end