%load all you images for  section segmentation in one folder. 
%This script has NO dapi input. no layer information needed. HOWEVER, the
%brainno is also automatically detected so you need to name your files with
%the 1st 4 numbers as brain no. If brain number is 43, name it 0043. 
%start. 
close all; 
clear all; 
i = dir('*.tif'); %lists all stk stacks in current directory
%Parameters you need to change
%%
Brainno = '005x'; %Brain number
Age = 'P21'; % age of animal
Genotype = 'SST-cre'; % this is also the treatment condition
Region = 'SSCtx'; %type the region (e.g.,SSCtx, MCtx, VCtx )
R1 = 1; %reporter channel; 
C1 = 2; %TF/quantifying channel; 


%%

for k = 1:numel(i) % loops through all images in the directory
    filename = i(k).name;
    Brainno = extractBefore(filename,5); 
    A = struct('Brainno', {Brainno}, 'Genotype', {Genotype}, 'Age', {Age}, 'Region', {Region});
    z=0;    
    while z == 0;       
        LRP_dach1(filename, R1, C1, A)
        questions = [{'Redo same images? Yes = 1, No = 2'}];
        default_answers = [{'2'}];
        answer = inputdlg(questions,'Redo',1,default_answers);
        a = str2double(answer(1)); 
        close all;
        if a ==2
            z = 1;
        else
            z = 0;
        end
    end

end

