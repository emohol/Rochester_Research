function [dp,c] = dprime(h,fa,hN,faN)

% 8/9/18
% adapted from DPRIME_SIMPLE by Karin Cox
% https://www.mathworks.com/matlabcentral/fileexchange/47711-dprime_simple-m?focused=3837057&tab=function
% 
%   inputs: 
%   hit = hit rate as float , or # of hits if hN is provided 
%   fA  = false alarm rate as float, or # of false alarms if faN is provided 
%   hN  = number of stimulus presentations
%   faN = number of lures, i.e. trials w/o a stimulus  
%
%   outputs: 
%   dp = d'
%   c = criterion c (negative values --> bias towards yes responses)

narginchk(2,4);

% check for values out of bounds
if or(or(h>1,h<0),or(fa>1,fa<0))
    if nargin == 4
        h  = h  / hN;
        fa = fa / faN;
    else
        error('input arguments must fall in the 0 to 1 range unless N is provided')
    end
end


if fa == 0 
    % implement standard correction
    fa = 1/(2*faN);
end
if h == 1
        % implement standard correction
    h = 1 - (1 / (2*hN)); 
end

if fa < 0.01 
    % implement standard correction
    fa = 0.01;
end

if h < 0.01
        % implement standard correction
    h = 0.01; 
end

%could implement special case for when h==0 d-prime = 0        
    
% d prime = z(h)-z(fA)
dp = norminv(h)-norminv(fa);

% c = -0.5*[z(h)+z(fA)]
c = -0.5*(norminv(h)+ norminv(fa));

end