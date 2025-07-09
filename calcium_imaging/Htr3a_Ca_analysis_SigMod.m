% Used for analysis of 2p in vivo calcium imaging data from HTR3a mice
% performing a head-fixed forlimb reaching task
%
% Requires: 
% "Data" structure created using Htr3a_CaImaging_analysis_part1
%
% This script creates a results_corr structure that has data about number
% of significantly modulated neurons and average baseline and movement Ca2+
% values per mouse
%
% Author: Richard Roth (rhroth@stanford.edu)
% Date: 2025

%% average trials first and then use fluorescence intensity at each frame as n for statistical test
% Data = Data_Reaching; Behavior = "Reaching start"; Modifier = "all"; 
% Data = Data_ReachingS; Behavior = "Reaching start"; Modifier = "Success";
% Data = Data_ReachingF; Behavior = "Reaching start"; Modifier = "Fail";
% Data = Data_Chewing; Behavior = "Chewing"; Modifier = "all"; 
% Data = Data_Grooming; Behavior = "Grooming"; Modifier = "all"; 
% Data = Data_ReachingFtoS; Behavior = "Reaching start"; Modifier = "Fail matched to S";


%behavior timepoint t is at frame number 41 
t = 41;

t_baseline = t-40:t-20;
t_post = t+1:t+21;

CurrentRegion_all = [];
for i = 1:size(Data,2)

CurrentRegion=[];
for Neuron = 1:Data(i).NumCells

CurrentCell_all = squeeze(Data(i).EventCaData_ZScore_base(:,Neuron,:));
CurrentCell_base = mean(squeeze(Data(i).EventCaData_ZScore_base(t_baseline,Neuron,:)),2);
CurrentCell_post = mean(squeeze(Data(i).EventCaData_ZScore_base(t_post,Neuron,:)),2);

[h,p] = Htr3a_permutationTest_f(CurrentCell_base', CurrentCell_post', 10000);

CurrentCell_postvbase = mean(CurrentCell_post)./mean(CurrentCell_base);

CurrentRegion(Neuron,1:5) = [h,p,CurrentCell_postvbase',mean(CurrentCell_base),mean(CurrentCell_post)];

end
CurrentRegion_all = [CurrentRegion_all ; CurrentRegion];
Data(i).Modulated = CurrentRegion;
end

%%
uniqueMice = unique([Data.Mouse]);

results = struct('Mouse', {}, 'TotalNumTrials', {}, 'TotalNumCells', {}, 'AllModulated', {}, 'TotalInc', {}, 'TotalDec', {}, 'TotalNS', {});
% Loop through each unique mouse
for i = 1:length(uniqueMice)
    mouse = uniqueMice(i);
    
    % Find indices where Mouse matches the current unique mouse
    idx = [Data.Mouse] == mouse;
    
    totalNumTrials = sum([Data(idx).NumTrials]); % Sum the NumTrials for the current mouse
    totalNumCells = sum([Data(idx).NumCells]); % Sum the NumCells for the current mouse
    allModulated = vertcat(Data(idx).Modulated); % Concatenate the Modulated matrices for the current mouse
    
    num_inc = sum(allModulated(:,1) == 1  & allModulated(:,5) > allModulated(:,4));
    num_dec = sum(allModulated(:,1) == 1  & allModulated(:,5) < allModulated(:,4));
    num_ns = sum(allModulated(:,1) == 0);
    num_all = num_inc+num_dec+num_ns;
    
    
    % Store the result
    results(i).Mouse = mouse;
    results(i).TotalNumTrials = totalNumTrials;
    results(i).TotalNumCells = totalNumCells;
    results(i).AllModulated = allModulated;   
    results(i).TotalInc = num_inc;
    results(i).TotalDec = num_dec;
    results(i).TotalNS = num_ns;   
    
    results(i).TotalIncPercent = num_inc./num_all;
    results(i).TotalDecPercent = num_dec./num_all;
    results(i).TotalNSPercent = num_ns./num_all;

    results(i).BaseAvg = mean(allModulated(:,4),1,'omitnan');
    results(i).PostAvg = mean(allModulated(:,5),1,'omitnan');
end

%% Calculate adjusted p-values and remake results
results_corr = results;
for i = 1:length(uniqueMice)
numComparisons = nnz(~isnan(results(i).AllModulated(:,4)));
pValues = results(i).AllModulated(:,2);
if numComparisons == 0
    continue
else
pValues_corr = min(pValues * numComparisons, 1);
end
results_corr(i).AllModulated(:,2) = pValues_corr;
results_corr(i).AllModulated(:,1) = pValues_corr < 0.05;

% recalculate number of modulated cells
    num_inc = sum(results_corr(i).AllModulated(:,1) == 1  & results_corr(i).AllModulated(:,5) > results_corr(i).AllModulated(:,4));
    num_dec = sum(results_corr(i).AllModulated(:,1) == 1  & results_corr(i).AllModulated(:,5) < results_corr(i).AllModulated(:,4));
    num_ns = sum(results_corr(i).AllModulated(:,1) == 0);
    num_all = num_inc+num_dec+num_ns;

    results_corr(i).TotalInc = num_inc;
    results_corr(i).TotalDec = num_dec;
    results_corr(i).TotalNS = num_ns;  

    results_corr(i).TotalIncPercent = num_inc./num_all;
    results_corr(i).TotalDecPercent = num_dec./num_all;
    results_corr(i).TotalNSPercent = num_ns./num_all;
end


