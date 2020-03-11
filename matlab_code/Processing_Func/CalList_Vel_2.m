function list = CalList_Vel()

% stimulus duration: stimOff-saccOn.
% this is because the stimulus duration
% is calculated from the moment the eye lands.

    
   
        %uservar - what the experimenter saved (specified variables to be
        %stored)
        list = eis_readData([], 'x');
        list = eis_readData(list, 'y');
               
        
        list = eis_readData(list, 'trigger', 'blink');
        list = eis_readData(list, 'trigger', 'notrack');
        list = eis_readData(list, 'trigger','saccade1');
        list = eis_readData(list, 'trigger','msaccade1');
%         list = eis_readData(list, 'trigger','drift1');
        list = eis_readData(list, 'trigger','em_event1');
        
        list = eis_readData(list,'photocell');
        list = eis_readData(list,'spf');
        list = eis_readData(list, 'uservar','resp');
        list = eis_readData(list, 'uservar','present');
        
        list = eis_readData(list, 'uservar','eccentricity');
        list = eis_readData(list, 'uservar','contrast');
        list = eis_readData(list, 'uservar','backgroundContrast');        
        list = eis_readData(list, 'uservar','backgroundImage');
        list = eis_readData(list, 'uservar','responseTime');
        list = eis_readData(list, 'uservar','presTime');
        list = eis_readData(list, 'uservar','spatialFreq');
        list = eis_readData(list, 'uservar','velocity_on');
        list = eis_readData(list, 'uservar','velocity_off');
        
        list = eis_readData(list, 'uservar','pixelAngle');
        list = eis_readData(list, 'uservar','xoffset');
        list = eis_readData(list, 'uservar','yoffset');
        list = eis_readData(list, 'uservar','fixTime');
        
        list = eis_readData(list, 'uservar','cuePrf');%when the peripheral cue starts to be shown
        list = eis_readData(list, 'uservar','fixOn');
        list = eis_readData(list, 'uservar','cueCtr');
        list = eis_readData(list, 'uservar','saccOn');
        list = eis_readData(list, 'uservar','saccOff');
        list = eis_readData(list, 'uservar','flashOn');
        list = eis_readData(list, 'uservar','rampOff');
        list = eis_readData(list, 'uservar','stimOff');
        list = eis_readData(list, 'uservar','quit');
          
        
        
        
%         list = eis_readData(list, 'uservar','SubjectName');
%         list = eis_readData(list, 'uservar','Debug');
%         list = eis_readData(list, 'uservar','Trial');
%         
%         list = eis_readData(list, 'uservar','StabCond');
%         list = eis_readData(list, 'uservar','RightEyeImage');
%         list = eis_readData(list, 'uservar','LeftEyeImage');
%         list = eis_readData(list, 'uservar','SFCond');
%         list = eis_readData(list, 'uservar','EyeCond');
%         list = eis_readData(list, 'uservar','Psudo');
%         
%         list = eis_readData(list, 'uservar','StimulusDuration');
%         list = eis_readData(list, 'uservar','FixationImage');
%         list = eis_readData(list, 'uservar','FixationSize');
%         list = eis_readData(list, 'uservar','FusionImage');
%         
%         
%         list = eis_readData(list, 'uservar','RecalFrequency');
%         list = eis_readData(list, 'uservar','RecalIncrement');
%         list = eis_readData(list, 'uservar','RecalOffsetX1');
%         list = eis_readData(list, 'uservar','RecalOffsetY1');
%         list = eis_readData(list, 'uservar','RecalOffsetX2');
%         list = eis_readData(list, 'uservar','RecalOffsetY2');
%         
%         list = eis_readData(list, 'uservar','RightOriginX');
%         list = eis_readData(list, 'uservar','RightOriginY');
%         list = eis_readData(list, 'uservar','LeftOriginX');
%         list = eis_readData(list, 'uservar','LeftOriginY');
%         
%         list = eis_readData(list, 'uservar','ViewDistance');
%         list = eis_readData(list, 'uservar','HorizontalSize');
%         list = eis_readData(list, 'uservar','VerticalSize');
%         list = eis_readData(list, 'uservar','HorizontalResolution');
%         list = eis_readData(list, 'uservar','VerticalResolution');
%         list = eis_readData(list, 'uservar','RefreshRate');
%         list = eis_readData(list, 'uservar','PixelAngle');
           
        
   
     