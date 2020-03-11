A = imread('NaturalImages/2.jpg');
I = rgb2gray(A);
I2 = im2double(I);
fh=1;
figure(fh)
imshow(I)

% A = imread('NaturalImages/1.bmp');
% I2 = im2double(A);
% fh=1;
% figure(fh)
% imshow(I2)

% F = fft2(A);
% imagesc(abs(fftshift(F)))

% FFT of image and matrix inversion
IM1 = fftshift( fft2(I2));
% IM1 = fft2(I2);

% magnitude is set to one
% IM2 = IM1./abs(IM1);
IM2 = IM1;

% builds X and Y of distance from centre
[X, Y]=meshgrid( ...
    [ round(size(IM2,2)/2):-1:1 ...
    1:round(size(IM2,2)/2) ],...
    [ round(size(IM2,1)/2):-1:1 ...
    1:round(size(IM2,1)/2) ]);

R = sqrt(X.^2 + Y.^2);
fh=fh+1;
figure(fh)
imshow(IM2./abs(IM2),[])
colorbar

fh=fh+1;
figure(fh)
imshow(R,[])
colorbar

% mu = 600;
sigma = 100000; %High sigma -> High frequency pass
gsize = size(R);
[X1,Y1] = ndgrid(1:gsize(1), 1:gsize(2));
yc = round(size(IM2,2)/2);
xc = round(size(IM2,1)/2);
exponent = ((X1-xc).^2 + (Y1-yc).^2)./(2*sigma);
% exponent = ((X-xc).^2 + (Y-yc).^2)./(2*sigma);
amplitude = 1 / (sigma * sqrt(2*pi));  
% val       = amplitude*(exp(-exponent));
val = exp(-exponent);
% inv_val = 1./val;

% IM3 = IM2;
IM3 = IM2.*val;
% IM3 = IM2 .* inv_val;
fh=fh+1;
figure(fh)
imshow(IM3./abs(IM3),[])
colorbar

% n = real( ifft2(ifftshift(IM2)) );
n = real( ifft2(ifftshift(IM3)) );
% nm = n ./ sqrt(sum( n(:).^2 )) ;
filteredIm = n;
fh=fh+1;
figure(fh)
imshow(filteredIm,[])
colorbar


%% Low pass to high pass 
gsize = size(R);
df = 1;
% for sigma = 1:df:100000
for sigma = 1:df:50
    [X1,Y1] = ndgrid(1:gsize(1), 1:gsize(2));
    yc = round(size(IM2,2)/2);
    xc = round(size(IM2,1)/2);
    exponent = ((X1-xc).^2 + (Y1-yc).^2)./(2*sigma);
    amplitude = 1 / (sigma * sqrt(2*pi));  
%     val       = amplitude*(exp(-exponent));
    val = exp(-exponent);
    IM3 = IM2.*val;
    n = real( ifft2(ifftshift(IM3)) );
    filteredIm = n;
    figure(5)
    imshow(filteredIm,[])
%     colorbar
end

df = 5;
for sigma = 50:df:500
    [X1,Y1] = ndgrid(1:gsize(1), 1:gsize(2));
    yc = round(size(IM2,2)/2);
    xc = round(size(IM2,1)/2);
    exponent = ((X1-xc).^2 + (Y1-yc).^2)./(2*sigma);
    amplitude = 1 / (sigma * sqrt(2*pi));  
%     val       = amplitude*(exp(-exponent));
    val = exp(-exponent);
    IM3 = IM2.*val;
    n = real( ifft2(ifftshift(IM3)) );
    filteredIm = n;
    figure(5)
    imshow(filteredIm,[])
%     colorbar
end

df = 50;
for sigma = 500:df:10000
    [X1,Y1] = ndgrid(1:gsize(1), 1:gsize(2));
    yc = round(size(IM2,2)/2);
    xc = round(size(IM2,1)/2);
    exponent = ((X1-xc).^2 + (Y1-yc).^2)./(2*sigma);
    amplitude = 1 / (sigma * sqrt(2*pi));  
%     val       = amplitude*(exp(-exponent));
    val = exp(-exponent);
    IM3 = IM2.*val;
    n = real( ifft2(ifftshift(IM3)) );
    filteredIm = n;
    figure(5)
    imshow(filteredIm,[])
%     colorbar
end

df = 500;
for sigma = 10000:df:100000
    [X1,Y1] = ndgrid(1:gsize(1), 1:gsize(2));
    yc = round(size(IM2,2)/2);
    xc = round(size(IM2,1)/2);
    exponent = ((X1-xc).^2 + (Y1-yc).^2)./(2*sigma);
    amplitude = 1 / (sigma * sqrt(2*pi));  
%     val       = amplitude*(exp(-exponent));
    val = exp(-exponent);
    IM3 = IM2.*val;
    n = real( ifft2(ifftshift(IM3)) );
    filteredIm = n;
    figure(5)
    imshow(filteredIm,[])
%     colorbar
end
fprintf('DONE! \n \n')
%% Bunch of other things that are probably unimportant right now

test = normpdf(R,mu,sigma);
% % IM3 = IM2.*test;
% IM3 = imfilter(IM2,test);

% IM3 = IM2;
% B = abs(IM2(:))>0.8;
% if (any(abs(IM3(:))>0.8))
%     IM3
%
% end

% F = 1./R;
% F = R.*test;

F = imfilter(R,test);
% F = imgaussfilt(R,'FilterDomain','frequency');
IM3 = IM2 .* F;
fh=fh+1;
figure(fh)
imshow(IM3,[])
colorbar

n = real( ifft2(ifftshift(IM3)) );
nm = n ./ sqrt(sum( n(:).^2 )) ;
filteredIm = n;
fh=fh+1;
figure(fh)
imshow(filteredIm,[])


%%
R_1 = zeros(500,1000);
ind = R<600;
for i=1:size(ind,1)
    for j=1:size(ind,2)
        if ind(i,j)
            R_1(i,j)= R(ind(i,j));
        end
    end
end
IM3 = IM2 .* R_1;
fh=fh+1;
figure(fh)
imshow(IM3,[])
colorbar

n = real( ifft2(ifftshift(IM3)) );
nm = n ./ sqrt(sum( n(:).^2 )) ;
filteredIm = n;
fh=fh+1;
figure(fh)
imshow(filteredIm,[])

%%
sigma = 10;
for mu=0:20:500
    test = normpdf(R,mu,sigma);
    % % IM3 = IM2.*test;
    % IM3 = imfilter(IM2,test);

    % IM3 = IM2;
    % B = abs(IM2(:))>0.8;
    % if (any(abs(IM3(:))>0.8))
    %     IM3
    %     
    % end

    % F = 1./R;
    % F = R.*test;
    F = imfilter(R,test);
    IM3 = IM2 .* F;
%     fh=fh+1;
%     figure(fh)
%     imshow(IM3,[])
%     colorbar
    
    n = real( ifft2(ifftshift(IM3)) );
    nm = n ./ sqrt(sum( n(:).^2 )) ;
    filteredIm = n;
%     fh=fh+1;
    figure(5)
    imshow(filteredIm,[])
end

%%
% load the images
% A = imread('NaturalImages/1.jpg');
 images    = cell(30,1);
%  cmap = gray(256);
 [I, cmap] = imread('NaturalImages/1.bmp');
 for i=1:17
     images{i} = imread(sprintf('NaturalImages/%d.bmp',i));
 end
 % create the video writer with 30 fps
 writerObj = VideoWriter('Velocity.avi');
 writerObj.FrameRate = 1;
   % open the video writer
   open(writerObj);
   % write the frames to the video
    for u=1:17    
       % convert the image to a frame
       frame = im2frame(images{u},cmap);
       writeVideo(writerObj, frame);
%        for v=1:30
%            writeVideo(writerObj, frame);
%        end
   end
   % close the writer object
   close(writerObj);
%    implay('Velocity.avi');

%%
function mat = gauss2d(mat, sigma, center)
gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));
mat = gaussC(R,C, sigma, center);
end

function val = gaussC(x, y, sigma, center)
xc = center(1);
yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(2*sigma);
val       = (exp(-exponent));
end
%%
% figure(4)
% pwelch(I2)
% colorbar 

% figure(2)
% plot(X,R)
% colorbar
% 
% figure(3)
% plot(Y,R)
% colorbar
