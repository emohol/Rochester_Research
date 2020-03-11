clear; clf

noise_path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Lum\noise';
directory = dir([noise_path filesep '*.bmp']);
pixelAngle = 1.06; %Vis dynamics


for i=1:length(directory)
    filePath = fullfile(noise_path,directory(i).name);
    noise = imread(filePath,'bmp'); 
    noise = double(noise(:,:,1));
    %noise = (noise - min(min(noise))) ./ ( max(max(noise)) - min(min(noise)) );
    noise = (noise - 0) ./ ( 255 - 0) ;
    
    imsize = size(noise,1); 
    
    xx = meshgrid(-imsize/2:imsize/2-1);
    yy = xx';
    r = sqrt (xx.^2 + yy.^2) * pixelAngle / 60;
     
     ECC = [0 4 8]; 
     for ecc = 1:length(ECC)
         mask =  r <= ECC(ecc)+1;
          subplot(length(ECC),1,ecc);
         plot(noise(mask)); hold on 
         
         if ecc == 1
             xlimts = xlim; 
         else
             xlim(xlimts); 
         end
      

     end
end