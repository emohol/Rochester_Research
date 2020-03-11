path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';
% subject = 'Nikunj';
% subject = 'A092';
subject = 'A036';
pathtodata = fullfile(path,subject);

% fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures';
sub_fig_path = fullfile(fig_path,subject);
if (~isfolder(sub_fig_path))
    mkdir(sub_fig_path)
end

read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

[valid,counter] = countingTrialsNK(ppt);

[noFiles,~] = size(ppt);