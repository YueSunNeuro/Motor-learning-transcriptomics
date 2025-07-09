% Used for analysis of 2p in vivo calcium imaging data from HTR3a mice
% performing a head-fixed forlimb reaching task
%
% This function is necessary for script Htr3a_CaImaging_analysis_part1 
% 
% Author: Richard Roth (rhroth@stanford.edu)
% Date: 2025
function [EventCaData, EventCaData_ZScore, EventCaData_ZScore_base, EventCaData_dFF, EventCaData_TrialAvg, EventCaData_ZScore_TrialAvg, EventCaData_ZScore_base_TrialAvg, EventCaData_dFF_TrialAvg, EventDLCData, EventDLCData_TrialAvg] = Htr3a_Ca_analysis_f(Event_filepath, Ca_filepath, DLC_filepath, Time_pre, Time_post, Behavior, Modifier)

%%Read calcium data (measured from FIJI ROIs)
Ca_data = readtable(Ca_filepath);
Ca_data = Ca_data(:,2:end);
Ca_data = Ca_data(:,2:2:end); %remove area numbers
Ca_data = table2array(Ca_data);

%%Read reach event data
Event_data = readtable(Event_filepath);
Event_data.Time = str2double(Event_data.Time);
Event_data.Frame = round(Event_data.Time .* Event_data.FPS);

%%Read DLC data
DLC_data = readtable(DLC_filepath);
DLC_data = DLC_data(:,2:end);
DLC_data = DLC_data(:,10:12); %These are the columns of dominent paw X Y likelyhood
DLC_data = table2array(DLC_data);

% Remove DLC data points with low DLC "likelyhood"
DLC_data_X = DLC_data(:,1);
DLC_data_Y = DLC_data(:,2);
DLC_data_X(DLC_data(:,3)<0.90) = NaN;
DLC_data_Y(DLC_data(:,3)<0.90) = NaN;
DLC_data(:,1) = DLC_data_X;
DLC_data(:,2) = DLC_data_Y;

DLC_data(:,1) = movmean(DLC_data(:,1),2);
DLC_data(:,2) = movmean(DLC_data(:,2),2);

% Calculate velocity from trajectory
dt = 1./100; % time interval between data points (100hz)
dx = diff(DLC_data(:,1)); % change in X coordinates between data points
dy = diff(DLC_data(:,2)); % change in Y coordinates between data points
DLC_speed = sqrt(dx.^2 + dy.^2) / dt; % velocity at each time point

start2p_row = Event_data.Behavior=="2P Start";
start2p_frame = Event_data.Frame(start2p_row,:);

end2p_row = Event_data.Behavior=="2p End";
end2p_frame = Event_data.Frame(end2p_row,:);

% Calculate which 2p frame corresponds to a given frame from the behavior
% camera. 2p imaged at ~20hz and camera at ~100hz

% First determine number of frames
numcCamframes = end2p_frame-start2p_frame+1;
num2pframes = size(Ca_data,1);

% subtract cam frames before 2p onset and multiply by frame ratio
% result is rounded up to neirest integer to match 2p frame
%frameCam=30278;
%frame2p = ceil((frameCam-start2p_frame+1)*num2pframes/numcCamframes);

Event_data.FrameIn2P = ceil((Event_data.Frame-start2p_frame+1)*num2pframes/numcCamframes);


%%

if exist("Modifier")
  CurrBeh_row = Event_data.Behavior==Behavior & Event_data.Modifier1==Modifier;
else
  CurrBeh_row = Event_data.Behavior==Behavior;
end

if Behavior=="Chewing"
  CurrBeh_row = Event_data.Behavior==Behavior & Event_data.Status=="START";
end

if Behavior=="Grooming"
  CurrBeh_row = Event_data.Behavior==Behavior & Event_data.Status=="START";
end


CurrBeh_2pframe = Event_data.FrameIn2P(CurrBeh_row,:);
CurrBeh_DLCframe = Event_data.Frame(CurrBeh_row,:);

CurrBeh_2pframe = CurrBeh_2pframe(CurrBeh_2pframe<6000-Time_post*20); %Remove events that are too close to end
CurrBeh_2pframe = CurrBeh_2pframe(CurrBeh_2pframe>Time_pre*20); %Remove events that are too close to start

CurrBeh_DLCframe = CurrBeh_DLCframe(CurrBeh_2pframe<6000-Time_post*20); %Remove events that are too close to end
CurrBeh_DLCframe = CurrBeh_DLCframe(CurrBeh_2pframe>Time_pre*20); %Remove events that are too close to start

Frames_total_Ca = Time_pre*20 + Time_post*20 +1;
Frames_total_DLC = Time_pre*100 + Time_post*100 +1;

EventCaData=NaN(Frames_total_Ca,size(Ca_data,2),size(CurrBeh_2pframe,1));
EventCaData_ZScore=NaN(Frames_total_Ca,size(Ca_data,2),size(CurrBeh_2pframe,1));
EventCaData_ZScore_base=NaN(Frames_total_Ca,size(Ca_data,2),size(CurrBeh_2pframe,1));
EventCaData_dFF=NaN(Frames_total_Ca,size(Ca_data,2),size(CurrBeh_2pframe,1));

EventDLCData=NaN(Frames_total_DLC,size(DLC_speed,2),size(CurrBeh_DLCframe,1));

for i=1:size(CurrBeh_2pframe,1)
current2pframe = CurrBeh_2pframe(i);
current_EventCaData = Ca_data(current2pframe-Time_pre*20:current2pframe+Time_post*20,:);
EventCaData(:,:,i) = current_EventCaData;
EventCaData_ZScore(:,:,i) = zscore(current_EventCaData,0,1);
EventCaData_ZScore_base(:,:,i) = (current_EventCaData-mean(current_EventCaData(1:(Time_pre*20),:)))./std(current_EventCaData(1:(Time_pre*20),:));
EventCaData_dFF(:,:,i) = (current_EventCaData-mean(current_EventCaData(1:(Time_pre*20),:)))./mean(current_EventCaData(1:(Time_pre*20),:));

currentDLCframe = CurrBeh_DLCframe(i);
EventDLCData(:,:,i)=DLC_speed(currentDLCframe-Time_pre*100:currentDLCframe+Time_post*100,:);
end

EventCaData_TrialAvg = mean(EventCaData,3);
EventCaData_ZScore_TrialAvg = mean(EventCaData_ZScore,3);
EventCaData_ZScore_base_TrialAvg = mean(EventCaData_ZScore_base,3);
EventCaData_dFF_TrialAvg = mean(EventCaData_dFF,3);
EventDLCData_TrialAvg = mean(EventDLCData,3,'omitnan');

end