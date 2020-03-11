function [thresh, par, varargout ] = psyfit(lvl, hits, varargin)
% PSYFIT() fits a psychometric function to the data.
%
% thresh = psyfit(lvl, hits)
% [thresh par] = psyfit(lvl, hits, trials)
% thresh = psyfit( ... , 'Chance', .25,'Lapses', 'Auto')
% thresh = psyfit( ... , 'Log', 'Flip')
% thresh = psyfit( ... , 'Title', title)
% thresh = psyfit( ... , 'Thresh', thresholdLevel)
% [thresh par threshB parB] = psyfit(lvl, hits, 'Boots', 2000)
%
% RETURN
%
% 'thresh' threshold value (abscissa at the inflection point).
%   Threshold extrapolation is by default disabled (function
%   will return NaN, see below how to force extrapolation).
%
% 'par' fitting parameters: muP sigmaP and lapses used to fit
% the function
%
% 'threshB' N-long array of bootstrapped thresholds.
%
% 'parB' N-long array of bootstrapped parameters.
%
% INPUT
%
% 'lvl' indicates the stimulus level for a particular trial
%   or group of trials.
%
% 'hits' indicates whether the response to a particular trial
%   was correct. For grouped trials it indicates how many
%   correct trials were scored at a particular stimulus level.
%
% 'trials' indicates if that given trial should be included
%   in the analysis (1) or not (0). When For grouped trials,
%   it indicates how many trials were run at that particular
%   stimulus level.
%
% 'lvl', 'hits' and 'trials' must be arrays of the same size and
%   can refer to single trials or multiple trials grouped
%   by testing level.
%
% 'Chance' sets the lower asymptote of the function
%   (the chance level of the experiment)
%
% 'Thresh' allows specifying the probability level at which
%   the threshold is computed. This is set by default to 75%.
%
% 'Lapses' sets the higher asymptote of the function
%   (lapsing rate). When 'Chance' is set to 0, Lapses
%   determines both higher and lower asymptote of the curve.
%
% 'Log' fits a log-normal function to the data. Log fits
%   are performed using 3 free parameters, unless Lapsing
%   rate is specified different than 0.
%
% 'Extra' enables threshold extrapolation.
%
% 'Boots' bootstrap N threshold values from the fitted
%   function (see WICHMANN&HILL 2001 II). The number of
%   iteration for bootstrapping must be indicated (2000
%   iterations generally provide a good estimate).
%
% 'Flip' flips the function upside-down to fit data
%   where performance decreases with increasing levels
%   of the testing variable (noise for example).
%
% 'PlotOff' hide the plot of psychometric function.
%
% written by Marco Boi
% last modified,
% 06 June 2012
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% sets up some default values for params

plotOn=1;
% this vector indicates which trials are picked
% for the analysis. By default, all trials are
% included.
trials = ones(size(hits));
distType = 'Normal';
lambdaP = 0;
gammaP = .5;
isFlipped = 0;
isLog = 0;
canExtrapolate = 1;
repeats = 0;
fname = [];
plotHandle = 1;
xl = 'Auto';
yl = 'Auto';
label = sprintf('%s: ML fit to the data', distType);
% check varargin for any non-default input
ii=1;
while ii < length(varargin)+1
    if sum( size( varargin{ii})) && ~ischar( varargin{ii})
        trials = cell2mat(varargin(ii));
        ii=ii+1;
    elseif isequal(char(varargin(ii)),'Chance')
        gammaP = cell2mat(varargin(ii+1));
        ii=ii+2;
    elseif isequal(char(varargin(ii)),'Lapses')
        if  iscellstr(varargin(ii+1))
            lambdaP = - isequal(...
                char(varargin(ii+1)),'Auto');
        else
            lambdaP = cell2mat(varargin(ii+1));
        end
        ii=ii+2;
    elseif isequal(char(varargin(ii)), 'DistType')
        distType = varargin{ii+1};
        ii = ii + 2;
    elseif isequal(char(varargin(ii)),'Thresh')
        thrP= cell2mat(varargin(ii+1));
        ii=ii+2;
    elseif isequal(char(varargin(ii)),'Title')
        label = char(varargin(ii+1));
        ii=ii+2;
    elseif isequal(char(varargin(ii)),'Save')
        fname = char(varargin(ii+1));
        ii=ii+2;
    elseif isequal(char(varargin(ii)),'Boots')
        repeats = cell2mat(varargin(ii+1));
        ii=ii+2;
    elseif isequal(char(varargin(ii)),'Flip')
        isFlipped = 1;
        ii=ii+1;
    elseif isequal(char(varargin(ii)),'Log')
        isLog = 1;
        ii=ii+1;
    elseif isequal(char(varargin(ii)),'Extra')
        canExtrapolate = 1;
        ii=ii+1;
    elseif isequal(char(varargin(ii)),'PlotOff')
        plotOn = 0;
        ii=ii+1;
    elseif isequal(char(varargin(ii)),'PlotHandle')
        plotHandle = varargin{ii+1};
        ii = ii + 2;
    elseif isequal(char(varargin(ii)),'Xlim')
        xl = varargin{ii+1};
        ii = ii + 2;
    elseif isequal(char(varargin(ii)),'Ylim')
        yl = varargin{ii+1};
        ii = ii + 2;
    else
        error('Sorry, i don''t think i know that argument...');
    end
end

% by default, threshold is computed at the inflection point
if ~ exist('thrP','var')
    thrP = gammaP + .5*(1-gammaP);
end

% this makes sure arrays are horizontal
if size(lvl,1)>size(lvl,2)
    lvl=lvl'; end
if size(hits,1)>size(hits,2)
    hits=hits';end
if size(trials,1)>size(trials,2)
    trials=trials';end

% sets up the actual arrays used for the calculation.
% For log-normal fits, we don't consider points at
% 0 level as they might create problems with the
% likelihood calculation.
if isLog
    x = double( unique( lvl(trials>0 & lvl>0) ) );
else
    x = double( unique( lvl(trials>0 ) ) );
end
k = zeros(size(x));
n = zeros(size(x));
for ii=1:size(x,2)
    k(ii) = sum(  hits ( (x(ii)==lvl) & trials )  );
    n(ii) = sum( trials( (x(ii)==lvl) ) );
end

% here one could specify values to constrain
% the parameter space.
switch distType
    case 'Normal'
        if isLog
            startPt(1) = mean( log( x ) );
            startPt(2)= std( log( x ) ) ;
            limits(:,1) = [-Inf +Inf];
            limits(:,2) = [0 +Inf];
        else
            startPt(1) = mean(x);
            startPt(2)= std(x);
            limits(:,1) = [-Inf +Inf];
%             limits(:,1) = [0 +Inf];
            limits(:,2) = [0 +Inf];
        end
    case 'Weibull'
        if isLog
            % the 0.632 quantile of the weibull is alpha
            startPt(1) = quantile( log( lvl ), 0.632 );
            startPt(2)= std( log( x ) ) ; % not sure about doing log(log(x))
            limits(:,1) = [0 +Inf];
            limits(:,2) = [0 +Inf];
        else
            startPt(1) = quantile(lvl, 0.632);
            startPt(2)= abs(regress(log(-log(1 - k./n+eps))', log(x)'));
            limits(:,1) = [0 +Inf];
            limits(:,2) = [0 +Inf];
        end
end

% here i can report the different models...
if lambdaP<0
    psyF = @(x, par) psyfun( x, par(1), par(2), gammaP,...
        par(3),isFlipped, isLog, distType);
    % starting point and range for lambdaP parameter
    limits(:,3) = [0 0.1];
    startPt(3)= .05;
else
    psyF = @(x, par) psyfun( x, par(1), par(2), gammaP,...
        lambdaP, isFlipped, isLog, distType);
end

% fit experiemntal data
[thresh , par] = findParams(k, canExtrapolate);

% probability at tested levels from fitted model
pi = psyF(x, par);

% thresholds and parameters bootstrapping
if repeats
    threshB = zeros(1,repeats);
    parB = zeros(length(par),repeats);
    for ii=1:repeats
%     parfor ii=1:repeats
        % bootstrapped hits
        kB = binornd( n , pi );
        [threshB(ii), parB(:,ii)] = findParams(kB, canExtrapolate);
        
%         if ~ mod(ii,10)
%             disp(['Bootstrapping...',num2str(ii)]);
%         end
    end
    varargout{ 1} = threshB;
    varargout{ 2} = parB;
end

% display warning message when threshold is zero
if isnan(thresh)
    msg = sprintf(['WARNING, the threshold value is not '...
        'comprised in the tested range.\n You can '...
        'obtain a threshold value by enabling '...
        'Extrapolation,\n but the estimate would not be accurate.']);
    disp(msg);
end

% goodness of fit as specified in WICHMANN&HILL 2001
g = gof(k,n,pi);
% chisq test as specified in WICHMANN&HILL 2001
chisq = chi2cdf( sum(  n.*((pi - k./n).^2) ./((pi+eps).*(1-pi+eps)) ) ,length(n)-1 );

% this is contained in a function so that can be re-used
% for bootstrapping the threshold.
    function [thresh , par] = findParams(k, canExtrapolate)

        % we make the likelihood negative so that it can be
        % minimized through fminsearch()
        lklF = @(par) -1 * likelihood(k, n, psyF(x, par), par, limits); %Zhetuo/NK: Could add offset to cost function here. 
        par = simplex( lklF, startPt);
        if lambdaP>=0
            par(3)=lambdaP;
        else
            lambdaP=par(3);
        end

        % inverse psychometric function
        thresh = invpsyfun( thrP, par(1), par(2), ...
            gammaP, lambdaP, isFlipped, isLog, distType);
        if ~canExtrapolate && (thresh > max(x) || thresh < min(x))
            thresh = NaN;
        end
    end


% do the drawing :)
if plotOn

    figure(plotHandle);

    % we take into account all trials
    x = double( unique( lvl(trials>0 ) ) );
    k = zeros(size(x));
    n = zeros(size(x));
    for ii=1:size(x,2)
        k(ii) = sum(  hits ( (x(ii)==lvl) & trials )  );
        n(ii) = sum( trials( (x(ii)==lvl) ) );
    end

    % for exponential fit, lower bound is set to 0
    if strcmp(xl, 'Auto') || strcmp(xl, 'auto')
        if isLog
            xLo = min(x)./sqrt( ( max( x)./min( x)));
            xHi = max(x).*sqrt( ( max( x)./min( x)));
        else
            xLo = min(x)-range(x)./2;
            xHi = max(x)+range(x)./2;
        end 
    else
        xLo = xl(1);
        xHi = xl(2);
    end
    xRg = xHi-xLo;

    if strcmp(yl, 'Auto') || strcmp(xl, 'auto')
        yLo = min([gammaP k./n]) -.05;
        yHi = 1.05;
    else
        yLo = yl(1);
        yHi = yl(2);
    end
    yRg = yHi - yLo;

    if ~ isnan(thresh)
        if isLog
            xi = exp( linspace( log( xLo),...
                log( xHi), 200));
        else
            xi = linspace(xLo, xHi, 200);
        end
        yi = psyF(xi,par);
        plot(xi,yi, 'LineWidth',2,'Color',[0 0.5 0.5] );
%         han=plot(xi,yi, 'LineWidth',2,'Color',[0 0.5 0.5] ); %NK
        hold on
    end
  %draws circles for displaying data points below
    for ii=1:length(x)
        plot( x(ii) , k(ii) / n(ii) ,...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(n(ii)/10));
        hold on
    end
    if isLog
        set( gca, 'xscale','log');
    end
 %SETS NUMBERS INSIDE DATA POINTS NK
    text(x , k ./ n, ...
        num2str(n'),'Color',[1 .0 .0],...
        'HorizontalAlignment','Center');

    axis( [xLo xHi yLo yHi] );

    % display text
    if isLog
        txtPosX = log( xHi-.3*xRg);
    else
        txtPosX = xHi-.3*xRg;
    end

    fontsize = 14;
    text( txtPosX , yLo +.4 * yRg,...
        strcat('threshold =',num2str(round(thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( txtPosX , yLo +.35 * yRg,...
        strcat('N = ',num2str(sum(n))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);

    if ~ isnan(thresh)

%         text(txtPosX , yLo +.30 * yRg,...
%             strcat('p({\chi^2}) = ',num2str(round(chisq*100)/100)),...
%             'FontSize',fontsize,'Color',[.3 .3 .5]);
%         text(txtPosX , yLo +.25 * yRg,...
%             strcat('p(D) = ',num2str(round(g*100)/100)),...
%             'FontSize',fontsize,'Color',[.3 .3 .5]);
        switch distType
            case 'Normal'
                text(txtPosX , yLo +.20 * yRg,...
                    strcat('{\mu} = ',num2str(round(par(1)*100)/100)),...
                    'FontSize',fontsize,'Color',[.3 .7 .5]);
                text(txtPosX , yLo +.15 * yRg,...
                    strcat('{\sigma} = ',num2str(round(par(2)*100)/100)),...
                    'FontSize',fontsize,'Color',[.3 .7 .5]);
            case 'Weibull'
                text(txtPosX , yLo +.20 * yRg,...
                    strcat('{\alpha} = ',num2str(round(par(1)*100)/100)),...
                    'FontSize',fontsize,'Color',[.3 .7 .5]);
                text(txtPosX , yLo +.15 * yRg,...
                    strcat('{\beta} = ',num2str(round(par(2)*100)/100)),...
                    'FontSize',fontsize,'Color',[.3 .7 .5]);
        end
        text(txtPosX , yLo +.10 * yRg,...
            strcat('{\lambda} = ',num2str(round(par(3)*100)/100)),...
            'FontSize',fontsize,'Color',[.3 .7 .5]);
        text(txtPosX , yLo +.05 * yRg,...
            strcat('{\gamma} = ',num2str(round(gammaP*100)/100)),...
            'FontSize',fontsize,'Color',[.3 .7 .5]);
        title(label,'FontSize',12);

        % Here we draw a line to mark the threshold
        line( [thresh thresh] , [thrP-lambdaP yLo], ...
            'color','r','LineStyle',':', 'LineWidth',2);
        line( [xLo thresh] , [thrP-lambdaP thrP-lambdaP],...
            'color','r','LineStyle',':', 'LineWidth',2);
    end
    hold off
    if ~isempty( fname)
        fname = [fname '.eps'];
        saveas( plotHandle, fname, 'epsc');
        saveas(plotHandle, [fname(1:end-4), '.fig'], 'fig');
    end

end
end