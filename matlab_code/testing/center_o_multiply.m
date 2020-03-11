x = linspace(0,2*pi,20); 
amp = 0.7;
% Incorrect way
figure(1)
imagesc(amp*(sin(x) + 1) /2); colormap gray; set(gca,'CLim',[0 1]); colorbar

% Correct way
figure(2)
imagesc(((amp*sin(x)) + 1) /2); colormap gray; set(gca,'CLim',[0 1]); colorbar