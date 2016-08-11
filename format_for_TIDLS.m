%Read a .txt file which has been output from python. Put this into a
%format which will be accepted by TIDLS codes. NOTE: the temperature
%should be exported in descending order to match this code. 
function format_for_TIDLS(data_filename,sample_param_filename,saveDir)
load(sample_param_filename); 
%Sort the structure to match what should have been loaded into the .txt
%file
sample_param = nestedSortStruct(python_inputs,'temperature'); 
[num_samples,dim] = size(sample_param); 
%Open the text file
fileID = fopen(data_filename); 
%Read line-by-line 
tline = fgetl(fileID);
line = 1; 
dataSave = cell(num_samples,1);  
sample_count = 1; 
new_line = 0; 
while ischar(tline)
    %If we have numbers, str2num should give us a 1D array
    this_line_num = str2num(tline); 
    if isempty(this_line_num) == 1
        %We either have the header or an empty line
        if line ~= 1
            %This is not the header, but signals that we are entering a new
            %line. Store the previous data and start a new matrix. 
            dataSave{sample_count} = data_now; 
            sample_count = sample_count+1; 
            %We only need to store the excess carrier density and lifetime
            new_line = 1; %signal that next time we make a new data_now vector
        end
    elseif isempty(this_line_num)==0
        if line == 2 || new_line == 1
            %We need to start the data_now storage
            data_now = [this_line_num(1,2) this_line_num(1,4)];
            new_line = 0; 
        else
            %Then we just add this line to our existing structure
            data_now(end+1,:) = [this_line_num(1,2) this_line_num(1,4)];
        end
    end
    tline = fgetl(fileID);
    %Increment the line counter
    line = line+1; 
end
%Save the very last sample 
dataSave{sample_count} = data_now; 
%close the file
fclose(fileID); 
info = sample_param;  
%Now we should have read the entire file. Save. 
save([saveDir '\Raw_data.mat'],'dataSave'); 
save([saveDir '\meas_info.mat'],'info'); 
