%Open two PL images, align them, and then find the lifetime ratio map. 
clear all; close all; 

filename = 'C:\Users\Mallory\Documents\Non-contact crucible\9-15-2015 experiment TR+Amanda\Processed data - all stages\Ingot bottom_Position 3_40LP_PLmaps.mat';
load(filename); 

% state1 = 'As-grown'; 
% state2 = 'STD'; 
% state3 = 'STD+TR'; 
% state4 = 'TR'; 
% state5 = 'TR+STD'; 

state1 = 'As-grown'; 
state2 = 'EXT'; 
state3 = 'EXT+TR'; 
state4 = 'TR'; 
state5 = 'TR+EXT'; 

% state1 = 'As-grown sister 1'; 
% state2 = 'As-grown sister 2'; 
% state3 = 'TR sister 1'; 
% state4 = 'TR sister 2'; 
% state5 = 'TR+EXT';
% state6 = 'TR+STD';

for i = 1:length(folders)
    if strcmp(folders{i},state1)==1
        index1 = i; 
    elseif strcmp(folders{i},state2)==1
        index2 = i; 
    elseif strcmp(folders{i},state3)==1
        index3 = i; 
    elseif strcmp(folders{i},state4)==1
        index4 = i; 
    elseif strcmp(folders{i},state5)==1
        index5 = i; 
    elseif strcmp(folders{i},state6)==1
    index6 = i; 
    end
end

sample_pairs = [index2,index3,0;index1,index4,index5];

%Compare each state to the base state (state1) 
% [state2_AG] = align_maps(store_cropped{index1},store_cropped{index2});
% figure; state2_AG_ratio = state2_AG./store_cropped{index1}; imagesc(state2_AG_ratio); axis('image');
% colorbar; title([state2 ' compared to As-grown']);
% [state3_AG] = align_maps(store_cropped{index1},store_cropped{index3});
% figure; state3_AG_ratio = state3_AG./store_cropped{index1}; imagesc(state3_AG_ratio); axis('image');
% colorbar; title([state3 ' compared to As-grown']);
% [state4_AG] = align_maps(store_cropped{index1},store_cropped{index4}); 
% figure; state4_AG_ratio = state4_AG./store_cropped{index1}; imagesc(state4_AG_ratio); axis('image');
% colorbar; title([state4 ' compared to As-grown']);
% [state5_AG] = align_maps(store_cropped{index1},store_cropped{index5}); 
% figure; state5_AG_ratio = state5_AG./store_cropped{index1}; imagesc(state5_AG_ratio); axis('image');
% colorbar; title([state5 ' compared to As-grown']);

% [state3_AG] = align_maps(store_cropped{index1},store_cropped{index3});
% figure; state3_AG_ratio = state3_AG./store_cropped{index1}; imagesc(state3_AG_ratio); axis('image');
% colorbar; title([state3 ' compared to As-grown']);
% [state4_AG] = align_maps(store_cropped{index2},store_cropped{index4});
% figure; state4_AG_ratio = state4_AG./store_cropped{index2}; imagesc(state4_AG_ratio); axis('image');
% colorbar; title([state4 ' compared to As-grown']);
% [state5_AG] = align_maps(store_cropped{index1},store_cropped{index5}); 
% figure; state5_AG_ratio = state5_AG./store_cropped{index1}; imagesc(state5_AG_ratio); axis('image');
% colorbar; title([state5 ' compared to As-grown']);
% [state6_AG] = align_maps(store_cropped{index2},store_cropped{index6}); 
% figure; state6_AG_ratio = state6_AG./store_cropped{index2}; imagesc(state6_AG_ratio); axis('image');
% colorbar; title([state6 ' compared to As-grown']);
% 
h=figure; 
imagesc(store_cropped{index1},[0 450]); 
axis('image');
colormap('gray');
colorbar; 
axis off; 
title([state1 ' 40LP']);
saveFile = [state1 '_40LP_scaled'];
hgsave(h,saveFile);
print(h,'-dpng',saveFile);

h=figure; 
imagesc(store_cropped{index2},[0 450]); 
axis('image');
colormap('gray');
colorbar; 
axis off; 
title([state2 ' 40LP']);
saveFile = [state2 '_40LP_scaled'];
hgsave(h,saveFile);
print(h,'-dpng',saveFile);

h=figure; 
imagesc(store_cropped{index5},[0 450]); 
axis('image');
colormap('gray');
colorbar; 
axis off; 
title([state5 ' 40LP']);
saveFile = [state5 '_40LP_scaled'];
hgsave(h,saveFile);
print(h,'-dpng',saveFile);


% figure; 
% imagesc(store_cropped{index5},[0 500]); 
% axis('image'); 
% colormap(gray); 
% 
% sample_fixed = imread('123-5_AG.tif');
% sample_moving = imread('123-5_TREXT.tif');
% cpselect(sample_moving,sample_fixed); 




