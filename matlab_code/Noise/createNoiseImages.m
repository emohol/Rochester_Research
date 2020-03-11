% generate random noise images to use
% need nnoise.m function from Utilities toolbox (see gitlab - aplabUR
% group)
rng(0); % for consistancy

noisedir = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code\main\noise';
% noisedir = 'noise';
if ~exist(noisedir, 'dir')
    mkdir(noisedir); 
end

nimages = 10;
imsize = 1920;
pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially

for ii = 1:nimages
    % save noise at full contrast
    noise = nnoise(imsize, pixelAngle);
    noise = noise / max(abs(noise(:))) / 2 + 1/2;
    
    fname = sprintf('noise_%i.bmp', ii-1);
    imwrite(repmat(noise, [1, 1, 3]), fullfile(noisedir, fname));%why the repmat here? to convert into RGB so that it's EYERIS readable
end