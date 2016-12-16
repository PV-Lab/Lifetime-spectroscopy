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

%Determine the wafer size based on the image
filename = 'C:\Users\Mallory\Documents\PERC mc-Si degradation\PERC LeTID May 25 2016\PERC LeTID047.tif';
raw = imread(filename); 
figure;
imshow(raw);
title('Raw data');
grayscale = rgb2gray(raw); 
figure;
imshow(grayscale); 
title('Grayscale data');
%Now threshold the image
thresh = im2bw(grayscale,0.2); 
figure;
imshow(thresh); 
title('Thresholded image');

%Number of items we're interested in
n = 7; 
for i = 1:n
    [thresh_cropped,x_crop,y_crop]=crop_PL(thresh); 
    %We want to know the number of pixels that have threshold == 0 (black)
    hold = find(thresh_cropped==0);
    [a,b]=size(hold);
    num_pixels(i) = a; 
end

%Now we need to convert the number of pixels to a size
conversion = 297.2/14039; %mm/pixel in one direction
area_pixel = conversion*conversion; %mm2 of one pixel since they are square

sample_face = num_pixels.*area_pixel; %mm2

