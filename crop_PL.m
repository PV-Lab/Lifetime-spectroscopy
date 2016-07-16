function [PLimage_now,x_crop,y_crop]=crop_PL(PLmap)
%crop the PL image
check = 0;
displayPL = figure; 
croppedPL = figure; 

while check == 0 
    PLimage_now = PLmap; 
    figure(displayPL);
    clf;
    imagesc(PLimage_now);
    colormap(gray);
    axis('image');

    %Get the coordinates for cropping
    disp('Choose the top left corner and then bottom right corner of the sample bounding box');
    [x_crop,y_crop] = ginput(2);
    [m,n] = size(PLimage_now); 
    %Get rid of highest y values we don't want
    PLimage_now(floor(y_crop(2)):m,:) = [];
    %Get rid of highest x values we don't want
    PLimage_now(:,floor(x_crop(2)):n) = [];
    %Get rid of lowest y values we don't want
    PLimage_now(1:floor(y_crop(1)),:) = [];
    %Get rid of lowest x values we don't want
    PLimage_now(:,1:floor(x_crop(1))) = [];
    figure(croppedPL);
    clf;
    imagesc(PLimage_now);
    axis('image');
    colormap(gray);

    promptcorrect='Is the cropping and apparent matching visually acceptable? Y/N';
    str = input(promptcorrect,'s');
    if isempty(str)
        str = 'Y';
    end

    if str == 'Y'
        check = 1;
    end
end