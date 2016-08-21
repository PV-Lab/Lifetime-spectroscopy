%This function writes data to .dat and .inf files which are compatible with
%QSSPL analyzer. 
function write_QSSPL_analyser(dirname,sample,T,matrix_to_write,Ctonum)
    filename_write = [dirname '\' sample '_' num2str(T) Ctonum 'averaged.Raw Data.dat'];
    %We need to write the first line
    first_line = '# Time (s)	Generation (V)	PC (V)	PL (V)';
    dlmwrite(filename_write,first_line,'delimiter',''); 
    %Now append the data
    dlmwrite(filename_write,matrix_to_write,'-append','delimiter','\t','precision','%-6.18e','newline','pc'); 
    %Make the .inf file
    copyfile([dirname '\' sample '_' num2str(T) Ctonum '1.inf'],[dirname '\' sample '_' num2str(T) Ctonum 'averaged.inf']);
end