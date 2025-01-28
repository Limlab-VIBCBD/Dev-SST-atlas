
function visplot2x1 (RG,  ch_seg)
r1 = imoverlay (RG, bwperim(ch_seg), [0 0 1]);

%get (hfig1, 'Position')
hfig1 = figure; 
set(hfig1, 'Position', [100 80 1200 850])

subplot_tight(1,2,1);
imshow(r1);
title ('reporter(GREEN) + Tdt(R) + circle '); 

subplot_tight(1,2,2);
imshow(ch_seg);
title ('GFP+ only '); 



