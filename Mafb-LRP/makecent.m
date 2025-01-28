%06/2016
%this function outputs the centroid of the bw image
%usage  is: centroidbw = makecent(yourbwimage); 

function centroidbw = makecent (ch_t)
stats1 = regionprops(ch_t, 'Centroid');
centroids1 = cat(1, stats1.Centroid);
cr1 = round (centroids1); 
x1 = cr1 (:, 1);
y1 = cr1 (:, 2); 
num_objects1 = length (stats1); 
indcent1 = sub2ind((size(ch_t)), y1, x1); 
centroidbw = false (size (ch_t)); 
centroidbw(indcent1) = true; 
end
