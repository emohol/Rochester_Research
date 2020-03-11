function f = RadialProfile(varargin)
% Evaluates the radial average of the matrix N*N
% starting from the center (N/2+1,N/2+1) and considering nsamples.
% Syntax: f = RadialProfile(img,nsamples)
%
% The optional parameter AngleStep specifies the angular increment
% (defaulty value 0.01)
%  Syntax: f = RadialProfile(img,nsamples,AngleStep)

img = varargin{1};
nsamples = varargin{2};

% if (nargin == 3) 
%     AngleStep = varargin(3);
% else
%     AngleStep = 0.01; 
% end
if nargin == 3
    if length(varargin{3}) > 1
        Angles = varargin{3};
    else
        Angles = 0:varargin{3}:2*pi; 
    end
else
    Angles = 0:.01:2*pi;
end


xcenter = floor(size(img,1)/2 + 1);
ccprof = zeros(1,nsamples);   % Performs the radial avaraging
ccdir = zeros(1,nsamples);
DirCounter = 0;
for theta = Angles
	for rad = 0:nsamples-1
		jj = round(rad*cos(theta));
		ii = round(rad*sin(theta));
		ccdir(rad+1) =  img(xcenter+ii,xcenter+jj);
	end
	DirCounter = DirCounter+1;
	ccprof = ccprof + ccdir;
end
ccprof = ccprof/DirCounter;
f = ccprof;