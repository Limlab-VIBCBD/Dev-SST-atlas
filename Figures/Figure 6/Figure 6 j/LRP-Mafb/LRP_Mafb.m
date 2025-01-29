%2023-10-04 this script segments the area with a reporter, 
% than quantify a second channel to decide on a cutoff intensity for 
% double positive cells
%R1 = reporter channel
%C1 = TF/quantifying channel

function LRP_Mafb (filename, R1, C1, A)

[~,imagename, ~] = fileparts(filename); 
%%
%structure filennames
s = struct('inputfile', {strcat(imagename, '.tif')},...
    'outputfile', {strcat(imagename, '_GFP-Mafb.csv')},...
    'outputraw', {strcat(imagename, '_RawInt.csv')},...
    'fig1', {strcat(imagename, '_fig1.png')},...
    'fig2', {strcat(imagename, '_fig2.png')},...
    'fig3', {strcat(imagename, '_fig3.png')},...
    'fig4', {strcat(imagename, '_fig4.png')},...
    'fig5', {strcat(imagename, '_fig5.png')}); 

%%
%getting ROI from reporter channel
[ch3_total2,  ~] = tdt_seg_20x(filename, R1, s);
%%
%%
%intensity calculation
[IntTF, BG, Ich1_s ] = cal_int(s, C1, ch3_total2);% col1 = TF intensity, col3 = row-index
[IntYFP, ~, Ich2_s ] = cal_int(s, R1, ch3_total2);% col1 = YFP intensity, col3 = row-index
%%
figure, 
hold on,
hplot = histfit(IntTF(:, 1), 10, 'kernel');

curve = hplot(2);
xC= get(curve, 'XData'); % Get x-coordinates of the fitted curve
yC = get(curve, 'YData'); % Get y-coordinates of the fitted curve

[maxY, idx] = max(yC); % Get the maximum y value and its index
maxX = xC(idx); % Get the corresponding x value of fit curve = mean MafB

% Calculate standard deviation from the curve
% Full Width at Half Maximum (FWHM) can be used to estimate sigma:
half_max = maxY / 2; % Half of the max y-value
indices = find(yC >= half_max); % Indices where y >= half_max
fwhm = xC(indices(end)) - xC(indices(1)); % Distance between these points
sigma_fit = fwhm / (2 * sqrt(2 * log(2))); % Convert FWHM to sigma - this is stdev of fitted curve

y=get(hplot,'YData');
x = get(hplot, 'Xdata'); 
datay = y{2, 1}; 
datax = x{2, 1}; 
clear x, clear y; 
diff = BG+ (0.57*sigma_fit) ; %2XBG
xline(diff,'-.y', 'LineWidth',2)
xline(maxX,'-.k', 'LineWidth',2)

text(max (datax)*0.6, max (datay)*0.8, sprintf('BG+0.57*STD', diff), 'FontSize', 12, 'Color', 'black');


hold off
saveas(gcf, s.fig2, 'png')
%clear IntTF

%%
Data1 = num2cell (IntTF);
%%
idx = IntTF(:,1)> diff;
Data1(:, 5) = num2cell (IntYFP (:, 1)); 
Data1(:, 4) = {'YFP-MafB-low'};%cell type annotation
Data1(idx,4) =  {'YFP-MafB-high'};%cell type

%%
T = cell2table(Data1);
T.Properties.VariableNames = ["MafB_int", 'index_col','index_row', 'Celltype', 'YFP_int']; 
T.Celltype = categorical(T.Celltype); 
T.Cellindex = (1:size(T,1)).'; 
T2 = T(T.Celltype == 'YFP-MafB-high', :); 
T3 = T(T.Celltype == 'YFP-MafB-low', :); 

%%
%display MafB
figure,imshowpair(imadjust(Ich1_s), ch3_total2, 'montage')
hold on, 
plot(T2.index_col, T2.index_row, 'ro', 'MarkerSize',5)
text(T2.index_col, T2.index_row*1.02, string(T2.Cellindex), 'Color','red','FontSize',10)
plot(T3.index_col, T3.index_row, 'bo', 'MarkerSize',5)
text(T3.index_col, T3.index_row*1.05, string(T3.Cellindex), 'Color','blue','FontSize',10)
title('MafB channel with YFP ROI')
hold off
saveas(gcf, s.fig3, 'png')
%%
%display YFP
figure,imshowpair(imadjust(Ich2_s), ch3_total2, 'montage')
hold on, 
plot(T2.index_col, T2.index_row, 'ro', 'MarkerSize',5)
text(T2.index_col, T2.index_row*1.02, string(T2.Cellindex), 'Color','red','FontSize',10)
plot(T3.index_col, T3.index_row, 'bo', 'MarkerSize',5)
text(T3.index_col, T3.index_row*1.02, string(T3.Cellindex), 'Color','blue','FontSize',10)
title('YFP channel with YFP ROI')
hold off
saveas(gcf, s.fig4, 'png')

%%
%display IntYPF and IntTF correlation:
figure, hold on
scatter (IntYFP (:, 1), IntTF (:, 1), 'b', 'filled')
%plot(polyfit(IntYFP (:, 1), IntTF (:, 1), 1), 'r'); % Line fit for clarity
yline(diff,'-.k', 'LineWidth',2)
xlabel('Int YFP)')
ylabel('Int Mafb')
title('Mafb vs YFP')
R = corrcoef(IntYFP (:, 1), IntTF (:, 1));
PearR = R (1, 2); 
text(max(IntYFP (:, 1)) * 0.7, max(IntTF (:, 1)) * 0.9, sprintf('PearsonR = %.2f', PearR), 'FontSize', 12, 'Color', 'black');
saveas(gcf, s.fig5, 'png')
hold off
%%
%R = corrcoef(IntYFP (:, 1), IntTF (:, 1));
G3 = groupcounts(T,["Celltype"],"IncludeEmptyGroups",true);
G3 = table2cell(G3);
G3(:, 3) = {A.Genotype}; 
G3(:, 4) = {A.Brainno}; 
G3(:, 5) = {A.Region}; 
G3(:, 6) = {imagename}; 

%%

Alldata = cell2dataset(G3, 'VarNames', {'Celltype', 'cellcount', 'Genotype', 'Brainno', 'Region','imagename'}); 
export (Alldata, 'File', s(1).outputfile,'Delimiter',',');
writetable (T, s(1).outputraw,'Delimiter',',');
end

