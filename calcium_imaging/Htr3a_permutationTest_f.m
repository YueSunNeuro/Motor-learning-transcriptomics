% Used for analysis of 2p in vivo calcium imaging data from HTR3a mice
% performing a head-fixed forlimb reaching task
%
% This function is necessary for script Htr3a_Ca_analysis_SigMod and performs 
% a two-tailed permutation test to determine whether the difference in means 
% between two input vectors (baseline and post movement initiation) is statistically significant.
% 
% Author: Richard Roth (rhroth@stanford.edu)
% Date: 2025
 
function [isSignificant, pValue] = Htr3a_permutationTest_f(vector1, vector2, numPermutations)
    % Combine the two vectors
    combined = [vector1, vector2];

    % Calculate the observed difference in means
    observedDiff = mean(vector1) - mean(vector2);

    % Initialize a counter for the number of permutations with a greater or equal difference
    numExtreme = 0;

    % Perform the permutation test
    for i = 1:numPermutations
        % Shuffle the combined data
        shuffled = combined(randperm(length(combined)));
        
        % Split the shuffled data into two groups
        permVector1 = shuffled(1:length(vector1));
        permVector2 = shuffled(length(vector1)+1:end);
        
        % Calculate the difference in means for the shuffled data
        permDiff = mean(permVector1) - mean(permVector2);
        
        % Check if the permuted difference is as extreme or more extreme than the observed difference
        if abs(permDiff) >= abs(observedDiff)
            numExtreme = numExtreme + 1;
        end
    end

    % Calculate the p-value
    pValue = numExtreme / numPermutations;

    % Determine if the result is significant
    isSignificant = pValue < 0.05;
end
