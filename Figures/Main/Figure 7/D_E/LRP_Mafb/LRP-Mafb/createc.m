%counting cells in centroid matrix. 
%cbw = the centroid image
%Ind = index file
%x = size of index file (row); 

function  = createc (cbw, Ind, x)
[~, t_count(1)] = bwlabel(cbw(1:Ind(1), :));  
    for i = 2:x
        [~, t_count(i)] = bwlabel(cbw(Ind(i-1):Ind(i), :)); 
    end
end