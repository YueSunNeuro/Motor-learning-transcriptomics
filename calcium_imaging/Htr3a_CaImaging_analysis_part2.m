% Used for analysis of 2p in vivo calcium imaging data from HTR3a mice
% performing a head-fixed forlimb reaching task
%
% Requires: 
% "Data" structure created using Htr3a_CaImaging_analysis_part1
%
% This script plots the "Data" structure created using Htr3a_CaImaging_analysis_part1
% 
% Author: Richard Roth (rhroth@stanford.edu)
% Date: 2025

%% Plot heatmaps 
% clearvars -except Data Behavior Modifier
load('Htr3aCalciumImagingFilenames.mat')
Time_pre = 2; %time in seconds before event
Time_post = 8; %time in seconds after event

DataSortBy = Data;
All_Neurons=[];
All_Neurons_SortBy=[];

% Use Baseline Z-Score
for i=1:size(Data,2)
    All_Neurons = [All_Neurons, Data(i).EventCaData_ZScore_base_TrialAvg];
    All_Neurons_SortBy = [All_Neurons_SortBy, DataSortBy(i).EventCaData_ZScore_base_TrialAvg];
end

All_Neurons_zscore_base = All_Neurons;
All_Neurons_SortBy_zscore_base = All_Neurons_SortBy;

figure

data_plot = All_Neurons_zscore_base';
[M,sortby] = max(All_Neurons_SortBy_zscore_base);
data_plot = [data_plot,sortby'];
[~, ses] = size(data_plot);
datasorted = sortrows(data_plot,ses);

datatoplot = datasorted(:,1:end-1);
datatoplot = removeNanRows(datatoplot);

[nr,nc] = size(datatoplot);
pcolor([datatoplot nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
clim([-0.5 2]);
colorbar;
title([Behavior Modifier "Calcium"])

ax = gca;

%% Plot average Ca2+ trace
x = 1:size(All_Neurons_zscore_base,1);
y = All_Neurons_zscore_base';

% Calculate the mean and 95% confidence interval of the data
y_mean = mean(y,1,'omitnan');
y_ci = tinv([0.025; 0.975], length(y)-1) * std(y,'omitnan') / sqrt(length(y));
y_lower = y_mean + y_ci(1);
y_upper = y_mean + y_ci(2);

% Create a shaded area using the fill function
figure
% subplot(3,1,2)
x_fill = [x, fliplr(x)];
y_fill = [y_upper, fliplr(y_lower)];
fill(x_fill, y_fill, [0.5, 0.7, 0.5], 'LineStyle', 'none');

hold on;

% Plot the mean as a line
plot(x, y_mean, 'LineWidth', 2, 'Color', 'g');
xlim([0 (Time_pre+Time_post)*20]);
ylim([-0.5 1]);
% ylim([-0.2 0.2]);
title([Behavior Modifier "Average Calcium"])

%% Plot average DLC forelimb trajectories 
All_DLC_trials=[];
for i=1:size(Htr3aCalciumImagingFilenames,1)
    All_DLC_trials = [All_DLC_trials, squeeze(Data(i).EventDLCData)];
end

x = 1:size(All_DLC_trials,1);
y = All_DLC_trials';


% Calculate the mean and 95% confidence interval of the data
y_mean = mean(y,1,'omitnan');
y_ci = tinv([0.025; 0.975], length(y)-1) * std(y,'omitnan') / sqrt(length(y));
y_lower = y_mean + y_ci(1);
y_upper = y_mean + y_ci(2);

% Create a shaded area using the fill function
figure
% subplot(3,1,3)
x_fill = [x, fliplr(x)];
y_fill = [y_upper, fliplr(y_lower)];
fill(x_fill, y_fill, [0.7, 0.7, 0.7], 'LineStyle', 'none');

hold on;

% Plot the mean as a line
plot(x, y_mean, 'LineWidth', 2, 'Color', 'k');
xlim([0 (Time_pre+Time_post)*100]);
ylim([0 3000]);
title([Behavior Modifier "Average Reach Velocity"])



%% Plot each mouse avg individually and global average
Time_pre = 2; %time in seconds before event
Time_post = 8; %time in seconds after event

uniqueMice = unique([Data.Mouse]);

data_mice = struct('Mouse', {}, 'TotalNumTrials', {}, 'TotalNumCells', {}, 'AllModulated', {});
% Loop through each unique mouse
for i = 1:length(uniqueMice)
    mouse = uniqueMice(i);
    
    % Find indices where Mouse matches the current unique mouse
    idx = [Data.Mouse] == mouse;
    
    totalNumTrials = sum([Data(idx).NumTrials]); % Sum the NumTrials for the current mouse
    totalNumCells = sum([Data(idx).NumCells]); % Sum the NumCells for the current mouse
    allNeurons = horzcat(Data(idx).EventCaData_ZScore_base_TrialAvg); % Concatenate all matrices for the current mouse
    allNeuronsAvg = mean(allNeurons,2);
    allTrials = squeeze(cat(3,Data(idx).EventDLCData)); % Concatenate all matrices for the current mouse
    allTrialsAvg = mean(allTrials,2,'omitnan');

    % Store the result
    data_mice(i).Mouse = mouse;
    data_mice(i).TotalNumTrials = totalNumTrials;
    data_mice(i).TotalNumCells = totalNumCells;
    data_mice(i).allNeurons = allNeurons;   
    data_mice(i).allNeuronsAvg = allNeuronsAvg;  
    data_mice(i).allTrials = allTrials;  
    data_mice(i).allTrialsAvg = allTrialsAvg;  
end

x = 1:size(data_mice(1).allNeurons,1);
y = [data_mice.allNeuronsAvg]';

figure
plot(x,y,'Color', [0.8, 0.8, 0.8])

% Calculate the mean and SEM of the data
y_mean = mean(y,1,'omitnan');
y_std = std(y,'omitnan');
y_sem = std(y,'omitnan') / sqrt(size(y,1));
y_lower = y_mean - y_sem;
y_upper = y_mean + y_sem;

% Create a shaded area using the fill function
% figure
hold on
% subplot(3,1,2)
x_fill = [x, fliplr(x)];
y_fill = [y_upper, fliplr(y_lower)];
fill(x_fill, y_fill, [0.5, 0.7, 0.5], 'LineStyle', 'none');

hold on;

% Plot the mean as a line
plot(x, y_mean, 'LineWidth', 2, 'Color', 'g');
xlim([0 (Time_pre+Time_post)*20]);
% ylim([-0.2 2]);
ylim([-1 2]);
% ylim([-0.2 0.2]);
title([Behavior Modifier "Average Calcium"])


% DLC trajectories

x = 1:size(data_mice(1).allTrials,1);
y = [data_mice.allTrialsAvg]';

figure
plot(x,y,'Color', [0.8, 0.8, 0.8])

% Calculate the mean and SEM of the data
y_mean = mean(y,1,'omitnan');
y_std = std(y,'omitnan');
y_sem = std(y,'omitnan') / sqrt(size(y,1));
y_lower = y_mean - y_sem;
y_upper = y_mean + y_sem;

% Create a shaded area using the fill function
% figure
hold on
x_fill = [x, fliplr(x)];
y_fill = [y_upper, fliplr(y_lower)];
fill(x_fill, y_fill, [0.5, 0.5, 0.5], 'LineStyle', 'none');
hold on;

% Plot the mean as a line
plot(x, y_mean, 'LineWidth', 2, 'Color', 'k');
xlim([0 (Time_pre+Time_post)*100]);
ylim([-500 6000]);
title([Behavior Modifier "Average Reach Velocity"])


