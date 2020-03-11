function [data, filenames] = readdata(pt, list, dosave, fname)

% LIST is the list of variables to import for 
% the specific experiment.
if ~exist('fname', 'var')
    fname = fullfile(pt, 'results.mat');%[pt,'/results.mat'];
end
if ~exist('dosave', 'var')
    dosave = false; 
end

% this converts a whole directory
data = eis_eisdir2mat(pt, list, fname);
filenames = dir(fullfile(pt, '*.eis'));

if dosave
    save(fname, 'data', 'filenames');
end

% if fname
%     load(fullfile(pt, 'results.mat'));
% end
