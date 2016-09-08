%This script is for obtaining the raw data from a Sinton spreadsheet,
%preparing it for calculation of lifetime from these values. 
function [unprocessed] = get_Sinton_raw(dirname)
    [fileList,fileListShort] = getAllFiles(dirname); 
    for file = 1:length(fileList)
        filename = fileList{file}; 
        %Grab the data
        data = xlsread(filename,'RawData','A2:C126'); 
        time{file,1} = data(:,1); 
        photovolt{file,1} = data(:,2); 
        refvolt{file,1} = data(:,3); 
        %Let's plot this raw data just for visual confirmation
        figure;
        plot(time{file,1},refvolt{file,1},'b-','LineWidth',2); 
        hold on;
        plot(time{file,1},photovolt{file,1},'r-','LineWidth',2); 
        xlabel('Time [s]','FontSize',20); 
        ylabel('Voltage [V]','FontSize',20); 
        title(fileListShort{file},'FontSize',20);
        legend('Reference','Photovoltage'); 
        %Let's also grab the relevant sample parameters
        thick{file,1} = xlsread(filename,'User','B6');
        res{file,1} = xlsread(filename,'User','C6');
        oc{file,1} = xlsread(filename,'User','E6');
        meas_res{file,1} = xlsread(filename,'Summary','Q2');
        doping{file,1} = xlsread(filename,'Summary','E2');
        [x,txt] = xlsread(filename,'User','D6'); 
        type{file,1} = txt{1}; 
        %Get the measurement parameters
        [x,txt] = xlsread(filename,'Summary','I2');
%         mode{file,1} = txt{1}; 
        mode{file,1} = '1/1';
        %Get the calibration parameters
        calib_A{file,1} = xlsread(filename,'Settings','C6');
        calib_B{file,1} = xlsread(filename,'Settings','C7');
        calib_offset{file,1} = xlsread(filename,'Settings','C8');
        calib_refcell{file,1} = xlsread(filename,'Settings','C5');
    end
    unprocessed = struct('filename',fileListShort,'thickness',thick,...
    'resistivity',res,'measured_resistivity',meas_res,...
    'optical_constant',oc,'doping',doping,'photovoltage',...
    photovolt,'ref_voltage',refvolt,'time',time,'A',calib_A,'B',...
    calib_B,'offset',calib_offset,'ref_cell_conversion',...
    calib_refcell,'type',type,'measurement_mode',mode);
    %Save the raw data in the same directory
    save([dirname '\voltage_data.mat'],'unprocessed','fileListShort'); 
end