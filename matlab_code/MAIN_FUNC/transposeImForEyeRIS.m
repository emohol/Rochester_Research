imgdir = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code\NaturalImages1';

Imgs = dir(fullfile(imgdir,'*.jpg'));


for j=1:length(Imgs)
    Img = imread(fullfile(imgdir,Imgs(j).name));  % Read image
    Im = double(transpose(Img(:,:,1)));
    
    Im1 = Im / max(abs(Im(:)));%range(0,1)
    
    fname = sprintf('trans_%i.bin',j);
    fname = fullfile(imgdir, fname);
    fid = fopen( fname,'w');
    fwrite( fid,Im1,'float');
    fclose( fid);
    
    fname = sprintf('trans_%i.bmp', j);
    imwrite(Im1, fullfile(imgdir,fname));%why the repmat here? to convert into RGB so that it's EYERIS readable
%     figure    
%     imshow(Im)
%     pause
end