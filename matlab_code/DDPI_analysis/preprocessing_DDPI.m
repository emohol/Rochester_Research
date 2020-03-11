% 
% This file is a modified version of Preprocess_Data()
% It takes a variable containing a list of the entries
% to implement.
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process eye movements data and create a list of trials.
% Only the period of stimulus presentation is examined.
% Eye movements are classified for all trials. No selection of trials is
% operated.
%
%  [ValidTrials] = PreprocessTrials(data,NumberOfTrials,MinSaccSpeed, MinSaccAmp,...
%                            MinMSaccSpeed, MinMSaccAmp, MaxMSaccAmp,RefractoryPeriod)


function [ValidTrials] = preprocessing_DDPI(data)
%     data,NumberOfTrials,MinSaccSpeed, MinSaccAmp,...
%     MinMSaccSpeed, MinMSaccAmp, MaxMSaccAmp,MinVelocity);

% THESE PARAMS SHOULD BE DEFAULT TO BE OVERWRITTEN ONLY
% WHEN A NEW SET OF PARAMS IS SPECIFIED
 
MinSaccSpeed = 180;
MinSaccAmp = 30;
MinMSaccSpeed = 180;
MinMSaccAmp = 3;
MaxMSaccAmp = 30;
MinVelocity = 120;
NumberOfTrials = length(data.x);

fprintf('\t Preprocessing trials\n');

ValidTrials = {};
ValidTrialsCounter = 0;

for TrialNumber = 1:NumberOfTrials
    DataValid = [];
    % Get the data relative to the trial

    X = (data.x{TrialNumber});
    Y = (data.y{TrialNumber});

    Blink = (data.triggers{TrialNumber}.blink);
    NoTrack = (data.triggers{TrialNumber}.notrack);
    DataValid = ones(1,length(data.x{TrialNumber}));
    Invalid = zeros(size(X));
    
    Trial = createTrial(TrialNumber,...
        X, Y, Blink, NoTrack, DataValid,Invalid,331);

    Trial = preprocessSignals(Trial, 'minvel',MinVelocity, 'noanalysis'); 
    %%%% 'noanalysis' keeps the saccades at the very onset of the trial
    
    
    % Find all valid saccades with speed greater than 3 deg/sec
    % and bigger than 30 arcmin
    Trial = findSaccades(Trial, 'minvel',MinSaccSpeed, 'minsa',MinSaccAmp);

    % Find all valid microsaccades with speed greater than 3 deg/sec
    % and amplitude included in 3 arcmin and 60 arcmin
    Trial = findMicrosaccades(Trial, 'minvel',MinMSaccSpeed,'minmsa', MinMSaccAmp,'maxmsa', MaxMSaccAmp);
    Trial = findDrifts(Trial);

%     if any(data.stream00{TrialNumber}.data)
%         Xpos = data.stream00{TrialNumber}.data; %data
%         Ypos = data.stream01{TrialNumber}.data;
%         XT = data.stream00{TrialNumber}.ts; %times
%         YT = data.stream01{TrialNumber}.ts;
%         D = double(XT(end))-length(X);
%         if XT(1)==XT(2)
%             XT = XT(2:end);
%             YT = YT(2:end);
%             Xpos = Xpos(2:end);
%             Ypos = Ypos(2:end);
%         end
%         % here add five ms at the beginning
%         XStab = eis_expandVector(Xpos, XT, length(X), 'replica');
%         YStab = eis_expandVector(Ypos, YT, length(X), 'replica');
%     else
%         XStab = [];
%         YStab = [];
%     end
    %%%%%%% JI: TO DO note hardcoded pixel angle
%     Trial.XStab = XStab(1:end).*1.3989; %arcmin
%     Trial.YStab = YStab(1:end).*1.3989;
    
    if isfield(data, 'user')
        fn=fieldnames(data.user{TrialNumber});
        for ff=1:size(fn,1)
            entry = char(fn(ff));
            Trial.(entry) = data.user{TrialNumber}.(entry);
        end
    end
    
    if ~isempty(Trial)
        ValidTrials(ValidTrialsCounter+1, 1) = {Trial};
        ValidTrialsCounter = ValidTrialsCounter + 1;
    end

end
