% Used for analysis of 2p in vivo calcium imaging data from HTR3a mice
% performing a head-fixed forlimb reaching task
%
% Requires:
% - Fluorescence traces of GCaMP signal from neurons extracted using ImageJ/FIJI
% - CSV files from analysis of head-fixed reachin task with Deep Lab Cut
% - xlsx files from event detection using Boris
% - overview.mat file with filenames of these files for each mouse/imaging session
%
% This script creates a Data structure that can be used with
% Htr3a_CaImaging_analysis_part2 for plotting
%
% Author: Richard Roth (rhroth@stanford.edu)
% Date: 2025

%%
clearvars Modifier Data

% Behavior = "Reaching start"; Modifier = "all"; 
% Behavior = "Grooming"; Modifier = "all"; 
% Behavior = "Reaching start"; Modifier = "Success";
% Behavior = "Reaching start"; Modifier = "Fail";
Behavior = "Chewing"; Modifier = "all"; 
% Behavior = "Reaching start"; Modifier = "Fail matched to S";

load('Htr3aCalciumImagingFilenames.mat')
for i=1:size(Htr3aCalciumImagingFilenames,1)

Mousename = Htr3aCalciumImagingFilenames.Mouse(i);
    
Ca_filename = Htr3aCalciumImagingFilenames.Ca_Filepath(i);
Ca_filepath = join(["C:\Users\Richard\Dropbox\SU_Ding Lab\Htr3a\CaImagingResultsCopy\" Mousename "\" Ca_filename],'');
DLC_filename = Htr3aCalciumImagingFilenames.DLC_Filepath(i);
DLC_filepath = join(["G:\My Drive\Htr3a\HTR3A-headfixed_DLC\" Mousename "\dlc files\" DLC_filename],'');
Event_filename = Htr3aCalciumImagingFilenames.Event_Filepath(i);
Event_filepath = join(["G:\My Drive\Htr3a\HTR3A-headfixed_DLC\" Mousename "\" Event_filename],'');   

Time_pre = 2; %time in seconds before event
Time_post = 8; %time in seconds after event

if exist("Modifier")
[EventCaData, EventCaData_ZScore, EventCaData_ZScore_base, EventCaData_dFF, EventCaData_TrialAvg, EventCaData_ZScore_TrialAvg, EventCaData_ZScore_base_TrialAvg, EventCaData_dFF_TrialAvg, EventDLCData, EventDLCData_TrialAvg] = Htr3a_Ca_analysis_f(Event_filepath, Ca_filepath, DLC_filepath, Time_pre, Time_post, Behavior, Modifier);
else
[EventCaData, EventCaData_ZScore, EventCaData_ZScore_base, EventCaData_dFF, EventCaData_TrialAvg, EventCaData_ZScore_TrialAvg, EventCaData_ZScore_base_TrialAvg, EventCaData_dFF_TrialAvg, EventDLCData, EventDLCData_TrialAvg] = Htr3a_Ca_analysis_f(Event_filepath, Ca_filepath, DLC_filepath, Time_pre, Time_post, Behavior);
end

Data(i).Mouse = Mousename;
    Session = convertStringsToChars(Event_filename);
Data(i).Session = Session(1:12);
Data(i).EventCaData = EventCaData;
Data(i).EventCaData_ZScore = EventCaData_ZScore;
Data(i).EventCaData_ZScore_base = EventCaData_ZScore_base;
Data(i).EventCaData_dFF = EventCaData_dFF;
Data(i).EventCaData_TrialAvg = EventCaData_TrialAvg;
Data(i).EventCaData_ZScore_TrialAvg = EventCaData_ZScore_TrialAvg;
Data(i).EventCaData_ZScore_base_TrialAvg = EventCaData_ZScore_base_TrialAvg;
Data(i).EventCaData_dFF_TrialAvg = EventCaData_dFF_TrialAvg;
Data(i).EventDLCData = EventDLCData;
Data(i).EventDLCData_TrialAvg = EventDLCData_TrialAvg;
Data(i).NumTrials = size(EventCaData,3);
Data(i).NumCells = size(EventCaData,2);

end