%segmenting reporter channel 10x images, pass the background subtracted
%image and the se file to the function. 

function [RCE_total2 , RCE_total] = tdt_seg_20x_auto(filename, RCE_ch, s)
%%
%define struct element
se(1) = strel('line', 10, 5);
se(2) = strel('line', 10, 40);
se(3) = strel('disk', 2, 8);
se(4) = strel('disk', 3, 8);
se(5) = strel ('disk', 1, 8);
se(6) = strel('disk', 10);
se(7) = strel('line', 3, 0);
se(8) = strel('line', 5, 0);
%%
Ich1 = imread(filename, RCE_ch); %channel for RCE
%Ich1_m = medfilt2 (Ich1);
level = graythresh(Ich1); 
bw_ch1 = imbinarize(Ich1,level/2);


%%
RCE_total = imerode(bw_ch1, [se(2), se(1)]); 
%figure, imshow(RCE_total); 
%%
%RCE_total = imfill (bw_ch1, 'holes'); 
RCE_total = bwareaopen (RCE_total, 10);

%%
%%
[RCE_total2, cXy] = makecentxy(RCE_total); 
RCE_total2 = imdilate(RCE_total2, se(6)); 
figure,imshowpair(imadjust(Ich1), RCE_total2, 'montage')
hold on, 
plot(cXy(:, 1), cXy(:, 2), 'ro', 'MarkerSize',10, 'LineWidth', 2)
hold off

saveas(gcf, s.fig1, 'png')

end

