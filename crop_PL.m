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