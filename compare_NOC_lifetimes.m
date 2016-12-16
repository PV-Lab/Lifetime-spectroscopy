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

clear all; close all; 

file_location = 'C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Processed data - all stages';
ingot_heights = {'Ingot bottom','Ingot top'};
wafer_positions = {'Position 1','Position 2','Position 3','Position 4','Position 5','Position 6'};
scales = [[0 500];[0 2500]]; 
displayPL = figure; 
croppedPL = figure; 
crop_locations = cell(length(ingot_heights),length(wafer_positions)); 
crop_file = 'C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Processed data - all stages\cropping_parameters.mat';
linear_tau = cell(length(ingot_heights),length(wafer_positions)); 
%Do this for each ingot height
for height = 1:length(ingot_heights)
    %Do this for each position
    for position = 1:length(wafer_positions)
        
        %We want a lifetime figure for each position
        lifetime = figure('units','normalized','outerposition',[0 0 1 1]);
        %We specifically want a SRH lifetime curve for each position
        SRHlifetime = figure('units','normalized','outerposition',[0 0 1 1]); 
        %We want a PLI figure for each position
        photolum = figure('units','normalized','outerposition',[0 0 1 1]); 
        %We want a PLI figure for each position
        SRH_linear = figure('units','normalized','outerposition',[0 0 1 1]); 
        
        file_location_now = [file_location '\' ingot_heights{height} '\' wafer_positions{position}];
        %Get the list of folders in that folder
        folders_hold = dir(file_location_now); 
        %We only care about the strings that describe the states
        folders = {};
        for entries = 3:length(folders_hold)
            folders{end+1} = folders_hold(entries).name; 
        end
        
        %Now go through the different folders and grab and plot the data
        crop_locations_hold = cell(length(folders)); 
        store_cropped = cell(length(folders)); 
        linear_tau_hold = cell(length(folders));
        for state = 1:length(folders)
            file_location_now = [file_location '\' ingot_heights{height} '\' wafer_positions{position} '\' folders{state}];
            %Get the lifetime data
            files_hold = dir(file_location_now); 
            %We only care about the files that are relevant
            files = {}; 
            for entries = 3:length(files_hold)
                files{end+1} = files_hold(entries).name; 
            end
            %assign the indices to lifetime versus PCPL
            for i = 1:length(files)
                if isempty(strfind(files{i},'PCPL')) == 0
                    %If not empty, this is the PCPL file
                    PCPL = files{i};
                elseif isempty(strfind(files{i},'tauBreakdown')) == 0
                    %If not empty, this is the injection dependent file
                    tau_all = files{i}; 
                end
            end
            
            %Load the data finally!
            load([file_location '\' ingot_heights{height} '\' wafer_positions{position} '\' folders{state} '\' PCPL]); 
            load([file_location '\' ingot_heights{height} '\' wafer_positions{position} '\' folders{state} '\' tau_all]); 
            
            %Plot the injection dependent data
            figure(lifetime); 
            loglog(deltanq,tauq.*1e6,'LineWidth',3); 
            hold all; 
            figure(SRHlifetime); 
            loglog(deltanq,tau_SRH.*1e6,'LineWidth',3); 
            hold all; 
            
            %Get the linearized SRH lifetime
            if isempty(strfind(files{1},'123'))==0
                N_dop = 4.1e15;
            elseif isempty(strfind(files{1},'124'))==0
                N_dop = 3.8e15;
            elseif isempty(strfind(files{1},'19'))==0 || isempty(strfind(files{1},'20'))==0
                N_dop = 2.8e15;
            else
                disp('Could not find doping level');
            end
            %Get sample parameters at specified temperature
            T = 300;
            type = 'n'; 
            [Efi,Efv,p0,n0,Eiv] = adv_Model_gen(T,N_dop,type); 
            %Normalized carrier density
            if type == 'p'
                X = (n0+deltanq)./(p0+deltanq);
            elseif type == 'n'
                X = (p0+deltanq)./(n0+deltanq);
            end
            figure(SRH_linear); 
            plot(X,tau_SRH.*1e6,'LineWidth',3); 
            hold all; 
            
            linear_tau_hold{state} = [X;tau_SRH];
            if LP(1) ~= 40
                disp([ingot_heights{height} ' ' wafer_positions{position} ' ' folders{state} ': Not the right laser power!']);
            end
            if false
                %Crop the image
                [PLimage_now,x_crop,y_crop]=crop_PL(PLmaps{1});            
                crop_locations_hold{state} = [x_crop,y_crop];
            else
                PLimage_now = PLmaps{1}; 
                load(crop_file); 
                xy_crop = crop_locations{height,position}; 
                xy_crop = xy_crop{state}; 
                x_crop = xy_crop(:,1);
                y_crop = xy_crop(:,2); 
                [m,n] = size(PLimage_now); 
                %Get rid of highest y values we don't want
                PLimage_now(floor(y_crop(2)):m,:) = [];
                %Get rid of highest x values we don't want
                PLimage_now(:,floor(x_crop(2)):n) = [];
                %Get rid of lowest y values we don't want
                PLimage_now(1:floor(y_crop(1)),:) = [];
                %Get rid of lowest x values we don't want
                PLimage_now(:,1:floor(x_crop(1))) = [];
                store_cropped{state} = PLimage_now; 
            end
            figure(photolum); 
            subplot(1,length(folders),state); 
            imagesc(PLimage_now,scales(height,:));
            axis('image');
            colormap('gray'); 
            title(folders{state},'FontSize',30); 
            axis off;       
        end
        if false
            crop_locations{height,position} = crop_locations_hold; 
        end
        linear_tau{height,position} = linear_tau_hold; 
        %Make final adjustments to the figures and save
        figure(lifetime);
        xlabel('Excess carrier density (cm^-^3)','FontSize',40);
        ylabel('Lifetime (\mus)','FontSize',40);
        legend(folders);
        title([ingot_heights{height} ': ' wafer_positions{position}],'FontSize',40); 
%         hgsave(lifetime,[ingot_heights{height} '_' wafer_positions{position} '_Lifetime']);
%         print(lifetime,'-dpng','-r0',[ingot_heights{height} '_' wafer_positions{position} '_Lifetime.png']); 
        
        figure(SRHlifetime);
        xlabel('Excess carrier density (cm^-^3)','FontSize',40);
        ylabel('SRH Lifetime (\mus)','FontSize',40);
        legend(folders);
        title([ingot_heights{height} ': ' wafer_positions{position}],'FontSize',40); 
%         hgsave(SRHlifetime,[ingot_heights{height} '_' wafer_positions{position} '_SRHLifetime']);
%         print(SRHlifetime,'-dpng','-r0',[ingot_heights{height} '_' wafer_positions{position} '_SRHLifetime.png']); 
        
%         figure(photolum); 
%         hgsave(photolum,[ingot_heights{height} '_' wafer_positions{position} '_PLI']);
%         print(photolum,'-dpng','-r0',[ingot_heights{height} '_' wafer_positions{position} '_PLI.png']); 
        
        figure(SRH_linear); 
        xlabel('Normalized excess carrier density [-]','FontSize',40);
        ylabel('Lifetime [\mus]','FontSize',40);
        legend(folders);
        title([ingot_heights{height} ': ' wafer_positions{position}],'FontSize',40); 
%         hgsave(SRH_linear,[ingot_heights{height} '_' wafer_positions{position} '_SRHLinear']);
%         print(SRH_linear,'-dpng','-r0',[ingot_heights{height} '_' wafer_positions{position} '_SRHLinear.png']); 
        
%         save([ingot_heights{height} '_' wafer_positions{position} '_40LP_PLmaps.mat'],'store_cropped','folders');
        
        
    end
end

% save('cropping_parameters.mat','crop_locations'); 
            

