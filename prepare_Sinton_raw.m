%Prepare raw data for processing in Python. Dirname is the directory where
%the .mat files are saved. T is a vector which contains the sample
%temperatures that correspond to each file in the order they have been
%read. 
function prepare_Sinton_raw(dirname,T)
%Load the dark voltage data 
load([dirname '\dark_voltages.mat']);
%Load the lifetime data
load([dirname '\voltage_data.mat']);
%We will do the same thing for each file
for i = 1:length(fileListShort)
    photovolt = unprocessed(i).photovoltage; 
    refvolt = unprocessed(i).ref_voltage; 
    time = unprocessed(i).time; 
    dark_voltage = Vd(i); 
    %Add a row to the end as specified
    time(end+1) = 0; 
    refvolt(end+1) = 0; 
    photovolt(end+1) = dark_voltage; 
    PL = zeros(size(time)); 
    filename = fileListShort{i}; 
    filename = filename(1:end-5); 
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],{'# Time(s)'},'Sheet1','A1'); 
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],{'Generation (V)'},'Sheet1','B1'); 
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],{'PC (V)'},'Sheet1','C1');
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],{'PL (V)'},'Sheet1','D1');
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],time,'Sheet1','A2'); 
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],refvolt,'Sheet1','B2');
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],photovolt,'Sheet1','C2'); 
    xlswrite([dirname '\' filename '.Raw Data.xlsx'],PL,'Sheet1','D2');
    %We also want to save the other parameters
    temperature{i,1} = T(i); 
    %Calculate the reflectance based on the optical constant. 
    R{i,1} = 100*(1-unprocessed(i).optical_constant); 
    FS{i,1} = 0.038/unprocessed(i).ref_cell_conversion/(1.6e-19); 
    thickness{i,1} = unprocessed(i).thickness; 
    doping{i,1} = unprocessed(i).doping; 
    type{i,1} = unprocessed(i).type; 
end
python_inputs = struct('filename',fileListShort,'R',R,'FS',FS,...
    'temp',temperature,'thickness',thickness,'doping',doping,'type',type); 
save([dirname '\python_inputs.mat'],'python_inputs'); 