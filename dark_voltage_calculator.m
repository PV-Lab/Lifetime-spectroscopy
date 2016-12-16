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

%Determine the dark voltage from a Sinton measurement
function [Vd] = dark_voltage_calculator(dirname)
[fileList,fileListShort] = getAllFiles(dirname); 
for i = 1:length(fileList)
    filename = fileList{i}; 
    %We need the calibration constants
    calibration = xlsread(filename,'Settings','C6:C8'); 
%     calibration = xlsread(filename,'User','L6:N6'); %for old spreadsheet
    A = calibration(1); 
    B = calibration(2); 
    offset = calibration(3); 

    %We also need the sheet resistance calculated from the spreadsheet
    sheet_rho = xlsread(filename,'Summary','M2'); 
%     sheet_rho = xlsread(filename,'Summary','P2'); %for old spreadsheet

    %Now we calculate the dark voltage
    VdV0 = (-B+sqrt((B^2)+(4*A*(1/sheet_rho))))/(2*A); 
    Vd(i) = VdV0+offset; 
    sheet_resistance(i) = sheet_rho; 
    calibration_parameters{i} = [A,B,offset]; 
end
%Save the result
save([dirname '\dark_voltages.mat'],'Vd','sheet_resistance','calibration_parameters'); 
end