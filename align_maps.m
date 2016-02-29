function [moving_aligned]=align_maps(fixed,moving)

figure;
imshowpair(fixed,moving,'blend'); 
title('Maps before alignment'); 

[optimizer,metric] = imregconfig('Multimodal');

moving_aligned = imregister(moving,fixed,'affine',optimizer,metric); 

figure;
imshowpair(fixed,moving_aligned,'blend');
title('Maps after alignment');
