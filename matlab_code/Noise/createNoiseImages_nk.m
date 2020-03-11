% generate random noise images to use
% need nnoise.m function from Utilities toolbox (see gitlab - aplabUR
% group)
rng(0); % for consistancy

noisedir = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code\main\noise_nk(pos)';
% noisedir = 'noise';
if ~exist(noisedir, 'dir')
    mkdir(noisedir); 
end

nimages = 1000;
imsize_x = 1080;
imsize_y = 1920;
pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially

for ii = 1:nimages
    % save noise at full contrast
%     noise = nnoise(imsize, pixelAngle);
    noise = nnoise_nk(imsize_x,imsize_y, pixelAngle);
    noise = noise / max(abs(noise(:))) / 2 + 1/2;%range(0,1)
%     noise = noise / max(abs(noise(:)));%range (-1,1)
%     noise1 = repmat(noise, [1, 1, 3]);
    noise1 = noise;
    fname = sprintf('noise_%i.bin', ii-1);
%     fname = sprintf('noise_%i.fpi', ii-1);
    fname = fullfile(noisedir, fname);
    fid = fopen( fname,'w');
    fwrite( fid,transpose(noise1),'float');
%     fwrite( fid,noise,'float');
    fclose( fid);
    fname = sprintf('noise_%i.bmp', ii-1);
    imwrite(transpose(noise1), fullfile(noisedir, fname));%why the repmat here? to convert into RGB so that it's EYERIS readable
end

%To test the orientation of images between Matlab and C++
% noisedir = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code\main\noise_nk(no_tran)';
% if ~exist(noisedir, 'dir')
%     mkdir(noisedir); 
% end
% fname = sprintf('noise_test.bin');
% test_im = zeros(1080,1920);
% test_im(540:1080, 960:1920) = 1;
% test_im(1:540, 1:960) = 1;
% test_im(1:100,1:100) = 0.5;
% test_im(1:100,1820:1920) = 0.5;
% fname = fullfile(noisedir, fname);
% fid = fopen( fname,'w');
% fwrite( fid,transpose(test_im),'float');
% fclose( fid);
% fname = sprintf('noise_test.bmp');
% imwrite(transpose(test_im), fullfile(noisedir, fname));