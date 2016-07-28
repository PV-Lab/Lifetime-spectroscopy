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