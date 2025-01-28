function xy_cent = getxy_cent (ch_t)
stats1 = regionprops(ch_t, 'Centroid');
centroids1 = cat(1, stats1.Centroid);
cr1 = round (centroids1); 
x1 = cr1 (:, 1);
y1 = cr1 (:, 2); 
xy_cent = [x1, y1];
end