%Read a .txt file which has been output from a measurement and store the
%RAW DATA. Put this into a format which will be accepted by other scripts.
function [data_now,time,PC,PL,gen] = format_QSSPL_rawdata(data_filename,savename)
%Open the text file
fileID = fopen(data_filename); 
%Read line-by-line 
tline = fgetl(fileID);
line = 1; 
sample_count = 1; 
new_line = 0; 
while ischar(tline)
    %If we have numbers, str2num should give us a 1D array
    this_line_num = str2num(tline); 
    %Check whether this is the header. If not, proceed. 
    if isempty(this_line_num)==0
        if line == 2 || new_line == 1
            %We need to start the data_now storage
            data_now = [this_line_num(1,1) this_line_num(1,2) this_line_num(1,3) this_line_num(1,4)];
            new_line = 0; 
        else
            %Then we just add this line to our existing structure
            data_now(end+1,:) = [this_line_num(1,1) this_line_num(1,2) this_line_num(1,3) this_line_num(1,4)];
        end
    end
    tline = fgetl(fileID);
    %Increment the line counter
    line = line+1; 
end
%close the file
fclose(fileID); 
%Plot our results for the visual check
figure;
plot(data_now(:,1),data_now(:,2));
hold all;
plot(data_now(:,1),data_now(:,3));
hold all;
plot(data_now(:,1),data_now(:,4));
xlabel('Time [s]','FontSize',20); 
ylabel('Voltage [V]','FontSize',20);
legend('Generation','PC','PL'); 
time = data_now(:,1); 
PC = data_now(:,3); 
gen = data_now(:,2);
PL = data_now(:,4); 
dataSave{1} = data_now; 
%Now we should have read the entire file. Save. 
save(savename,'dataSave','PC','time','PL','gen');  