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

