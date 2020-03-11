function [L_FRs, LN_FRs] = RGC_LN_Model_Bin()
	%% Linear-nonlinear (LN) model for Retinal Ganglion Cells 

	sf = 3;			% spatial frequency
	swPix = 1366;	% screen width in pixels
	swMm = 600;		% screen width in mm
	shPix = 768;	% screen height in pixels
	shMm = 335;		% screen height in mm
	sDist = 1620;%1190;	% distance from screen to subject's eye in mm
	gwPix = 720;	% width of gabor patch in pixels; half of the screen height

	orientation = 0;
	phase = 0;
	wlPix = swPix ./ ( atand(swMm/2/sDist) * 2 * sf );
	contrast = 1;
	bgnLuminance = 128;
	img = Gabor( wlPix, orientation, phase, swPix, shPix, 'grating' )' * contrast;

	tRamp = 1500; % ms
    tPlateau = 1000;    % ms
	T = tRamp + tPlateau;	% ms
	Gt = ones(1,T);
	if(withRamp)
		Gt(1:tRamp) = (0:tRamp-1)/(tRamp-1);
	end

	% spatial receptive field
	rfWDeg = 5;%10;
	rfHDeg = 5;%10;
	rfWPix = ceil( tand(rfWDeg/2) * 2 * sDist / swMm * swPix );
	rfHPix = ceil( tand(rfHDeg/2) * 2 * sDist / shMm * shPix );
	precision = rfWDeg / rfWPix;
	rfLoc = round( ( [size(img)] - [rfWPix, rfHPix] ) / 2 ) - 1;	% bottom left point [x,y]
    rfImg = zeros(rfHPix, rfWPix, T);	% initialize visual input on rf
   
    
    vt = load( 'CharlesFV.mat', 'Trials' );
    Trials = vt.Trials;
    
	L_FRs = ones(size(Trials,2),T) * NaN;	% response from the  linear stage
	LN_FRs = ones(size(Trials,2),T) * NaN;	% firing rate from the whole LN model
	for( iTrial = 1:size(Trials,2) )
		% connec all drifts
		x = zeros( 1, sum([Trials(iTrial).drifts.duration]) );
		y = zeros(size(x));
		index = 0;
		for( i = 1 : size( Trials(iTrial).drifts.start, 2 ) )						
			idx = find( Trials(iTrial).drifts.start(i) == Trials(iTrial).blinks.start + Trials(iTrial).blinks.duration );	% whether drift(i) starts during a blink
			if( ~isempty(idx) )
				st = Trials(iTrial).drifts.start(i) + 250;
				dur = Trials(iTrial).drifts.duration(i) - 250;
				if( dur < 1 )
					continue;
				else
					Trials(iTrial).drifts.start(i) = st;
					Trials(iTrial).drifts.duration(i) = dur;
				end
			end
			idx = find( Trials(iTrial).drifts.start(i) + Trials(iTrial).drifts.duration(i) == Trials(iTrial).blinks.start );  % whether drift(i) finishes during a blink
			if( ~isempty(idx) )
				dur = Trials(iTrial).drifts.duration(i) - 50;
				if( dur < 1 )
					continue;
				else
					Trials(iTrial).drifts.duration(i) = dur;
				end
			end
			x( (1 : Trials(iTrial).drifts.duration(i)) + index ) = Trials(iTrial).x.position( (0 : Trials(iTrial).drifts.duration(i)-1) + Trials(iTrial).drifts.start(i) );
			y( (1 : Trials(iTrial).drifts.duration(i)) + index ) = Trials(iTrial).y.position( (0 : Trials(iTrial).drifts.duration(i)-1) + Trials(iTrial).drifts.start(i) );
			if( index > 0 )
				x( (1 : Trials(iTrial).drifts.duration(i)) + index ) = x( (1 : Trials(iTrial).drifts.duration(i)) + index ) - x(index+1) + x(index);
				y( (1 : Trials(iTrial).drifts.duration(i)) + index ) = y( (1 : Trials(iTrial).drifts.duration(i)) + index ) - y(index+1) + y(index);
			end
			index = index + Trials(iTrial).drifts.duration(i);
		end
		x( index+1 : end ) = [];
		y( index+1 : end ) = [];
		if( size(x,2) < T ) continue; end
		x = x(1:T);
		y = y(1:T);

		% rotate by -orientation
		xx = x * cosd(orientation) + y * sind(orientation);
		yy = y * cosd(orientation) - x * sind(orientation);

		% center eye trace
		x = x - mean(x);
		y = y - mean(y);

		% convert to pixels
		x = round( sDist * tand(x/60) / swMm * swPix );
		y = round( sDist * tand(y/60) / shMm * shPix );

		%% reconstruct visual input (moving under drifts)
		for( t = 1 : T )
			rfImg(:,:,t) = ( img( (1:rfWPix) + rfLoc(1) + x(t), (1:rfHPix) + rfLoc(2) + y(t) )' * Gt(t) + bgnLuminance );
		end

		% linear stage
		[fRate, fRate_C, fRate_S] = SpatioLinearModel( rfImg, precision, round(rfWPix/2), round(rfHPix/2), 'm', 25 );	% fRate: response from spatial RF; fRate_C: response from spatial RF center; fRate_S: response from spatial RF surround; fRate = fRate_C - fRate_S
        L_FRs(iTrial,:) = TemporalLinearModel( fRate, 'm', 'on' );		% L_FRs: response from the whole linear stage
		
		% Non-linear state: rectification
		base = 10;	% baseline activity
		LN_FRs(iTrial,:) = max( L_FRs(iTrial,:) + base, 0 );	% LN_FRs: firing rate from the whole LN model

	end
	L_FRs( isnan(L_FRs(:,1)), : ) = [];
	LN_FRs( isnan(LN_FRs(:,1)), : ) = [];
end


function [FR, FR_C, FR_S] = SpatioLinearModel( input, precision, xCenter, yCenter, cellType, eccentricity, isInterp )
	%% From Croner & Kaplan, Vision Research, 1995
    %% input:       2nd dimension as horizontal (x), 1st dimension as vertical (y); 3rd dimension as time
    %  precision:   degrees per pixel in input 
    %  xCenter:     horizontal location of the receptive field center on the input image, which is along the 2nd dimension of the input matrix
    %  yCenter:     vertical location of the receptive field cente on the input image, which is along the 1st dimension of the input matrx

    if( nargin() < 6 || eccentricity < 0 || eccentricity > 40 )
        index = 1;
        isInterp = false;   % interpolation for eccentricity infeasible in this case
    else
        index = 0;
    end
    if( nargin() < 7 )
        isInterp = false;   % no interpolation for eccentricity by default
    end

    if( lower(cellType) == 'p' )    % P cell
        ectRange = [ 0, 40; 0, 5; 5, 10; 10, 20; 20, 30; 30, 40 ];  % eccentricity ranges; degrees
        r_c = [ 0.05;   0.03;   0.05;   0.07;   0.09;   0.15 ]; % center radius, which gives 1/e of the peak sensitivity;   degrees
        K_c = [ 106.3;  325.2;  114.7;  77.8;   57.2;   18.6 ]; % center peak sensitivity;     spikes / ( s * %contrast * degree^2 )
        r_s = [ 0.42;   0.18;   0.43;   0.54;   0.73;   0.65 ]; % surround radius, which gives 1/e of the peak sensitivity;   degrees
        K_s = [ 1.1;    4.4;    0.7;    0.6;    0.8;    1.1 ];  % surround peak sensitivity;     spikes / ( s * %contrast * degree^2 )
    elseif( lower(cellType) == 'm' )    % M cell
        ectRange = [ 0, 40; 0, 10; 10, 20; 20, 30 ];  % eccentricity ranges; degrees
        r_c = [ 0.17;   0.10;   0.18;   0.23 ]; % center radius, which gives 1/e of the peak sensitivity;   degrees
        K_c = [ 84.7;   148.0;  115.0;  63.8 ]; % center peak sensitivity;     spikes / ( s * %contrast * degree^2 )
        r_s = [ 0.80;   0.72;   1.19;   0.58 ]; % surround radius, which gives 1/e of the peak sensitivity;   degrees
        K_s = [ 1.2;    1.1;    2.0;    1.6 ];  % surround peak sensitivity;     spikes / ( s * %contrast * degree^2 )
    else
        FR = [];
        FR_C = [];
        FR_S = [];
        return;
    end

    if( ~isInterp )
        if( index ~= 1 )
            index = find( ectRange(1:end,1) <= eccentricity & eccentricity <= ectRange(1:end,2), 1, 'first' );
            if( isempty(index) )
                return;
            else
                index = index(1);
            end
        end
        r_c = r_c(index);
        K_c = K_c(index);
        r_s = r_s(index);
        K_s = K_s(index);
    end

    x = ( (1 : size(input,2)) - xCenter ) * precision;
    y = ( (1 : size(input,1)) - yCenter ) * precision;
    [x,y] = meshgrid(x,y);
    x = x(:)';
    y = y(:)';
    input = reshape( input, [], size(input,3) );
    FR_C = K_c * exp( -(x.^2+y.^2) / r_c^2 ) * input * precision^2;
    FR_S = - K_s * exp( -(x.^2+y.^2) / r_s^2 ) * input * precision^2;
    FR = FR_C + FR_S;
end


function [FR, tRF ] = TemporalLinearModel( input, cellType, specifier )
	%% From Benardete & Kaplan, JPhy, 1999 (P cells); and Benardete & Kaplan, VisNeuro, 1999 (M cells)
	%  input:			firing rate resulted from spatial receptive field
	%  cellType:		P cell or M cell
	%  specifier:		ON cell or OFF cell when cellType is M; center or surround when cellType is P

   	FR = [];
   	tRF = [];
	
	if( lower(cellType) == 'p' )    % P cell
		if( strcmpi( specifier, 'center' ) )
            A = 13.98;     % spikes / (s * unit contrast)
            D = 4.0/1000;	% s, delay of transmission from the optic chiasm to the LGN
            N_L = 27;
            Tau_L = 2.05/1000;  % s
            H_S = 0.68;
            N_H = 1;
            Tau_H = 21.25/1000;	% s
        elseif( strcmpi( specifier, 'surround' ) )
        	A = 10.42;
        	D = 4.0/1000;
        	N_L = 37;
        	Tau_L = 1.54/1000;
        	H_S = 0.61;
        	N_H = 1;
        	Tau_H = 24.41/1000;
        else
        	return;
        end

    elseif( lower(cellType) == 'm' )    % M cell
    	if( strcmpi( specifier, 'on' ) )
    		A = 566.92;
            D = 2.2/1000;   % s
            N_L = 30.30;
            Tau_L = 1.41/1000;  % s
            H_S = 1;%0.98;        % H_S = 1 makes the summation of the kernel be zero
            Tau_0 = 54.6/1000;  % s
            C_half = 0.056;
        elseif( strcmpi( specifier, 'off' ) )
    		A = 550.14;
            D = 2.31/1000;   % s
            N_L = 22.60;
            Tau_L = 1.98/1000;  % s
            H_S = 1;%0.93;        % H_S = 1 makes the summation of the kernel be zero
            Tau_0 = 153.34/1000;  % s
            C_half = 0.051;
    	else
    		return;
    	end
    	N_H = 1;
    	c = 0.4;	% contrast
    	Tau_H = Tau_0 / ( 1 + (c/C_half)^2 );	% s
    
    else
    	return;
    end

	sRate = 1000;   % sampling rate of 1000 Hz
    nSamples = 1000;%round( .350 * sRate );   % number of samples for 350 ms; even number
    t = (1:nSamples) / sRate;    % 350 ms
    f = (0:nSamples/2) / (nSamples/sRate);
    w = 2*pi * f( 1 : nSamples/2+1 );
    K = A * exp( -i*w*D ) .* ( 1 - H_S ./ (1 + i*w*Tau_H) ).^N_H .* ( 1 ./ (1 + i*w*Tau_L) ).^N_L;
    K(end+1:end*2-2) = conj(K(end-1:-1:2));
    tRF = real(ifft(K));

    if(isempty(input)) return; end

    FR = conv( input, tRF );
    FR = FR(1:size(input,2));

end


function patch = Gabor( waveLength, orientation, pahse, width, height, window, sigma )
	%% patch = Gabor( waveLength, orientation, pahse, width, height, window, sigma )
    %  waveLength:      length of each cycle in pixels
    %  orientation:		counterclockwise; vertical gabor at 0; degrees
	%  pahse:			degrees
	%  width:			width of the patch in pixels (2nd dimension / x)
	%  height:			height of the patch in pixels (1st dimension / y)
	%  window:			'gaussian' (gabor), 'circle' (a round patch) or 'grating' (grating)
	%  sigma:			sigma of the Gaussian filter for 'gaussian' window, or of the Gaussian edge for 'circle' window (set to <=0 for a hard dented edge)
	%  patch:			1st dimension: vertical(y); 2nd dimension: horizontal(x)
	if( nargin() < 6 )
		window = 'grating';
	elseif( nargin() < 7 )
		sigma = 1;
	end
	
	[ x, y ] = meshgrid( (1:width) - (1+width)/2.0, (1:height) - (1+height)/2.0 );
	X = x.*cosd(orientation) + y.*sind(orientation);
	Y = y.*cosd(orientation) - x.*sind(orientation);
	frequency = 1/waveLength;
	if( strcmpi( window, 'gaussian' ) )
		patch = cos( 2 * pi * frequency .* X + pahse/180*pi ) .* exp( -0.5 * ((x/width*height).^2+y.^2) / sigma^2 );
	else
		patch = cos( 2 * pi * frequency .* X + pahse/180*pi );
	end
	if( strcmpi( window, 'circle' ) )
		mask = ones(height, width);
		if( sigma <= 0 )
			mask( ((x/width).^2+(y/height).^2) > 1/4 ) = 0;
		else
			index = (x/width*height).^2 + y.^2 >= (height/2-2*sigma)^2;			% make an edge of 2*sigma in vertical and 2*sigma/height*width in horizontal
			mask(index) = normpdf( sqrt( (x(index)/width*height).^2 + y(index).^2 ), height/2-2*sigma, sigma ) / normpdf(0,0,sigma);
		end
		patch = patch .* mask;
	end
	patch = patch / max(patch(:));
end