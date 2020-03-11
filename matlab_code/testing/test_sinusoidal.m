clf
count=1;
levels = [];
delta = 0.01;%.01
subplot(2,1,1)
for amplitude = [0.0:delta:0.5]   
%     amplitude = 0.02;
    x = sin(linspace(pi/2,3/2*pi,50));
    y = x*amplitude;
    y1 = ((y/2) + 0.5);
    y2 = round(y1*255);
    plot(y2);
    levels(1,count) = length(unique(y2));
    levels(2,count) = amplitude;
    levels(3,count) = median(y2);
    levels(4,count) = max(y2);
    levels(5,count) = min(y2);
    hold on
    count= count+1;
end
axis tight 
subplot(2,1,2)
plot([levels(2,:)],[levels(4,:)],'d-'); hold on
plot([levels(2,:)],[levels(5,:)],'d-'); hold on
xlabel('Signal Amp')
ylabel('Gray levels')  

if any(any(diff(levels(4:5,:)') == 0))
    title('bad step')
end