function [Int, m_SG, Ich1] = cal_int(s, C1, BW)
[M1, num_m] = bwlabel(BW); 
stats = regionprops(M1, 'PixelIdxList', 'Centroid'); 

% read images and background subtract
Ich1 = imread(s(1).inputfile, C1); %channel for TF 

for i=1:num_m 
   Int(i, 1)= mean(Ich1(stats(i).PixelIdxList));
   centroids1 = cat(1, stats(i).Centroid);
   Int(i, 2) = centroids1(1, 1); 
   Int(i, 3) = centroids1(1, 2); 
   
end
Int = round (Int); 
m_SG = mode (mean(Ich1));%backgound 
end

