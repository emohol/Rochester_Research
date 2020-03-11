function [span, mn_speed, mn_cur, varx, vary] = getDriftChar_nk(x, y, smoothing, cutseg, cutseg_front, maxSpeed)
% compute drift span
mx = mean(x(1+cutseg_front:end-cutseg));
my = mean(y(1+cutseg_front:end-cutseg));
span = max(sqrt((x - mx).^2 + (y - my).^2));

% s-golay filtering
smx = sgfilt(x, 3, smoothing, 0);
smy = sgfilt(y, 3, smoothing, 0);
smx1 = sgfilt(x, 3, smoothing, 1);
smy1 = sgfilt(y, 3, smoothing, 1);

% compute drift speed
sp = 1000*sqrt(smx1(1+cutseg_front:end-cutseg).^2 + smy1(1+cutseg_front:end-cutseg).^2);

if any(sp > maxSpeed)
    mn_speed = nan;
else
    mn_speed = mean(sp);
end


% compute drift curvature
tmpCur = abs(cur(smx(1+cutseg_front:end-cutseg), smy(1+cutseg_front:end-cutseg)));
mn_cur = mean(tmpCur);

varx = var(x(1+cutseg_front:end-cutseg));
vary = var(y(1+cutseg_front:end-cutseg));
end