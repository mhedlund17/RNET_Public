% Clean the workspace and close all figures
clear;
close all;

protocols = {'test'};
first_experiment_count = 0;
last_experiment_count = 34;
basePath = '/Users/sina/Documents/Stanford/Misc/D-Lab/mfp_analysis/data/';
mice = {'m1','m2', 'm3', 'm4', 'm5', 'm6'};

% Initialize photometry variables
numfibers = 2;
ITIlongthreshold = 19; %Threshold for fast ITI
ITIshortthreshold = 13; %Threshold for slow ITI

% Initialize combined variables
combinedPsth_struct = struct();
combinedRTToneUnsuccessful_struct = struct();
combinedRTToneSuccessful_struct = struct();
combinedRTUnsuccessful_struct = struct();
combinedRTSuccessful_struct = struct();
combinedPsthSuccess_struct = struct();
combinedPsthWithoutSecond_struct = struct();
combinedITISuccessful_struct = struct();
combinedITIUnsuccessful_struct = struct();
numtotaltrials = struct();
numsuccessfultrials = struct();
percentsuccessfultrials = struct();
combinednumtotaltrials = struct();
combinednumsuccessfultrials = struct();
combinedpercentsuccessfultrials = struct();
combinedPsthMatrix = [];
combinedRTToneUnsuccessful = [];
combinedRTToneSuccessful = [];
combinedRTUnsuccessful = [];
combinedRTSuccessful = [];
combinedPsthSuccessMatrix = [];
combinedPsthWithoutSecondMatrix = [];
combinedITISuccessful = [];
combinedITIUnsuccessful = [];

for protocolIdx = 1:length(protocols)
    protocol = protocols{protocolIdx};
    protocolPath = fullfile(basePath, protocol);

    for mouseIdx = 1:length(mice)
        mousePath = fullfile(protocolPath, mice{mouseIdx});
        addpath(mousePath);

        experimentNumbers = first_experiment_count:last_experiment_count;
        experiments = cell(1, length(experimentNumbers)); 
        for i = 1:length(experimentNumbers)
            experiments{i} = sprintf('experiment_%03d', experimentNumbers(i));  % Format the experiment number
        end

        combinedPsth_struct.(mice{mouseIdx}) = [];
        combinedRTToneSuccessful_struct.(mice{mouseIdx}) = [];
        combinedRTToneUnsuccessful_struct.(mice{mouseIdx}) = [];
        combinedPsthSuccess_struct.(mice{mouseIdx}) = [];
        combinedPsthWithoutSecond_struct.(mice{mouseIdx}) = [];
        combinedRTSuccessful_struct.(mice{mouseIdx}) = [];
        combinedRTUnsuccessful_struct.(mice{mouseIdx}) = [];
        combinednumtotaltrials.(mice{mouseIdx}) = [];
        combinednumsuccessfultrials.(mice{mouseIdx}) = [];
        combinedpercentsuccessfultrials.(mice{mouseIdx}) = [];
        combinedITISuccessful_struct.(mice{mouseIdx}) = [];
        combinedITIUnsuccessful_struct.(mice{mouseIdx}) = [];
    
        for experimentIdx = 1:length(experiments)
            experimentPath = fullfile(mousePath, experiments{experimentIdx});
            [psthMatrix, ITISuccessful, ITIUnsuccessful, RTToneSuccessful, RTToneUnsuccessful, RTSuccessful, RTUnsuccessful, psthSuccessMatrix, psthWithoutSecondMatrix, psthTime] = processDataForExperiment_zscored(experimentPath);
    
            % Combine/aggregate the outputs from each experiment
            combinedPsth_struct.(mice{mouseIdx}) = [combinedPsth_struct.(mice{mouseIdx}); psthMatrix];
            combinedRTToneSuccessful_struct.(mice{mouseIdx}) = [combinedRTToneSuccessful_struct.(mice{mouseIdx}), RTToneSuccessful];
            combinedRTToneUnsuccessful_struct.(mice{mouseIdx}) = [combinedRTToneUnsuccessful_struct.(mice{mouseIdx}), RTToneUnsuccessful];
            combinedITISuccessful_struct.(mice{mouseIdx}) = [combinedITISuccessful_struct.(mice{mouseIdx}), ITISuccessful];
            combinedITIUnsuccessful_struct.(mice{mouseIdx}) = [combinedITIUnsuccessful_struct.(mice{mouseIdx}), ITIUnsuccessful];
            combinedPsthSuccess_struct.(mice{mouseIdx}) = [combinedPsthSuccess_struct.(mice{mouseIdx}); psthSuccessMatrix];
            combinedPsthWithoutSecond_struct.(mice{mouseIdx}) = [combinedPsthWithoutSecond_struct.(mice{mouseIdx}); psthWithoutSecondMatrix];
            combinedRTSuccessful_struct.(mice{mouseIdx}) = [combinedRTSuccessful_struct.(mice{mouseIdx}), RTSuccessful];
            combinedRTUnsuccessful_struct.(mice{mouseIdx}) = [combinedRTUnsuccessful_struct.(mice{mouseIdx}), RTUnsuccessful];
            numtotaltrials.(mice{mouseIdx}).(experiments{experimentIdx}) = length(RTSuccessful);
            numsuccessfultrials.(mice{mouseIdx}).(experiments{experimentIdx}) = sum(~isnan(RTSuccessful));
            percentsuccessfultrials.(mice{mouseIdx}).(experiments{experimentIdx}) = 100 * (sum(~isnan(RTSuccessful))/length(RTSuccessful));
        end
    
         combinednumtotaltrials.(mice{mouseIdx}) = length(combinedRTSuccessful_struct.(mice{mouseIdx})) -1; %adjusting for removing the first trial
         combinednumsuccessfultrials.(mice{mouseIdx}) = sum(~isnan(combinedRTSuccessful_struct.(mice{mouseIdx})));
         combinedpercentsuccessfultrials.(mice{mouseIdx}) = 100 * (sum(~isnan(combinedRTSuccessful_struct.(mice{mouseIdx})))/length(combinedRTSuccessful_struct.(mice{mouseIdx})));
    end

    % Loop over each mouse and concatenate data to the combined matrices
    for p = 1:length(mice)
        mouseName = mice{p}; 
        
        % Extract data for the current mouse using dynamic field reference
        mouseData_combinedPsthMatrix = combinedPsth_struct.(mouseName);
        mouseData_combinedRTToneUnsuccessful = combinedRTToneUnsuccessful_struct.(mouseName);
        mouseData_combinedRTToneSuccessful = combinedRTToneSuccessful_struct.(mouseName);
        mouseData_combinedRTUnsuccessful = combinedRTUnsuccessful_struct.(mouseName);
        mouseData_combinedRTSuccessful = combinedRTSuccessful_struct.(mouseName);
        mouseData_combinedPsthSuccessMatrix = combinedPsthSuccess_struct.(mouseName);
        mouseData_combinedPsthWithoutSecondMatrix = combinedPsthWithoutSecond_struct.(mouseName);
        mouseData_combinedITISuccessful = combinedITISuccessful_struct.(mouseName);
        mouseData_combinedITIUnsuccessful = combinedITIUnsuccessful_struct.(mouseName);
    
        % Concatenate mouse data to the combined matrices along the first dimension
        combinedPsthMatrix = cat(1, combinedPsthMatrix, mouseData_combinedPsthMatrix);
        combinedRTToneUnsuccessful = [combinedRTToneUnsuccessful, mouseData_combinedRTToneUnsuccessful];
        combinedRTToneSuccessful = [combinedRTToneSuccessful, mouseData_combinedRTToneSuccessful];
        combinedRTUnsuccessful = [combinedRTUnsuccessful, mouseData_combinedRTUnsuccessful];
        combinedRTSuccessful = [combinedRTSuccessful, mouseData_combinedRTSuccessful];
        combinedPsthSuccessMatrix = cat(1, combinedPsthSuccessMatrix, mouseData_combinedPsthSuccessMatrix);
        combinedPsthWithoutSecondMatrix = cat(1, combinedPsthWithoutSecondMatrix, mouseData_combinedPsthWithoutSecondMatrix);
        combinedITISuccessful = [combinedITISuccessful, mouseData_combinedITISuccessful];
        combinedITIUnsuccessful = [combinedITIUnsuccessful, mouseData_combinedITIUnsuccessful];
    
    end
end

% check for habituation 
% Initialize variables to store RTs across all mice and experiments
allRTs_firstHalf = [];
allRTs_secondHalf = [];

% Loop over each mouse and each experiment to aggregate RTs
for mouseIdx = 1:length(mice)
    mouseName = mice{mouseIdx};
    experimentNumbers = 0:last_experiment_count;
    experiments = cell(1, length(experimentNumbers)); 
    for i = 1:length(experimentNumbers)
        experiments{i} = sprintf('experiment_%03d', experimentNumbers(i));  % Format the experiment number
    end
    
    for experimentIdx = 1:length(experiments)
        experimentName = experiments{experimentIdx};
        
        % Extract RTs for the current day
        RTToneSuccessful = combinedRTToneSuccessful_struct.(mouseName); % No need to use experimentName as we are aggregating for all experiments
        
        % Divide the trials into two halves
        numTrials = length(RTToneSuccessful);
        firstHalfRTs = RTToneSuccessful(1:floor(numTrials/2));
        secondHalfRTs = RTToneSuccessful(floor(numTrials/2)+1:end);
        
        % Append to the aggregated lists
        allRTs_firstHalf = [allRTs_firstHalf, firstHalfRTs];
        allRTs_secondHalf = [allRTs_secondHalf, secondHalfRTs];
    end
end

% Calculate mean and standard deviation for the first and second halves
meanFirstHalf = mean(allRTs_firstHalf);
meanSecondHalf = mean(allRTs_secondHalf);
stdFirstHalf = std(allRTs_firstHalf);
stdSecondHalf = std(allRTs_secondHalf);

% Perform a paired t-test across all aggregated data
[~, pval] = ttest2(allRTs_firstHalf, allRTs_secondHalf);

% Display aggregated results
disp('Aggregated Results Across All Trials:');
disp(['Mean First Half RT: ', num2str(meanFirstHalf)]);
disp(['Mean Second Half RT: ', num2str(meanSecondHalf)]);
disp(['P-value: ', num2str(pval)]);

% Plot the aggregated results
figure;
bar(1, meanFirstHalf, 'FaceColor', 'b'); hold on;
bar(2, meanSecondHalf, 'FaceColor', 'r');
errorbar(1, meanFirstHalf, stdFirstHalf, 'k', 'LineWidth', 2);
errorbar(2, meanSecondHalf, stdSecondHalf, 'k', 'LineWidth', 2);
xticks([1 2]);
xticklabels({'First Half', 'Second Half'});
ylabel('Mean RT (s)');
title('Comparison of Reaction Times between First and Second Half of Trials Across All Days');
text(1.5, max([meanFirstHalf + stdFirstHalf, meanSecondHalf + stdSecondHalf]) + 0.01, sprintf('p = %.3f', pval), 'HorizontalAlignment', 'center', 'FontSize', 12);
hold off;
box off;

% First, check if there are any successful trials for plotting.
if isempty(combinedPsthSuccessMatrix) || size(combinedPsthSuccessMatrix, 1) == 0
    fprintf('No successful trials to plot.\n');
else
    % Plot PSTHs for successful trials only
    figure(3);
    for p = 1:numfibers
        subplot(numfibers, 1, p);
        plot(psthTime, combinedPsthSuccessMatrix(:,:,p)');
        title(['Successful Trials PSTH - Fiber ', num2str(p), '- N = ', num2str(size(combinedRTToneSuccessful, 2)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
    end
    
    % Plot mean PSTH with error bars for successful trials only
    figure(4);
    for p = 1:numfibers
        subplot(numfibers, 1, p);
        plot(psthTime, mean(combinedPsthSuccessMatrix(:,:,p)));
        hold on;
        errSuccess = std(combinedPsthSuccessMatrix(:,:,p)) / sqrt(size(combinedPsthSuccessMatrix, 1));
        fill([psthTime, fliplr(psthTime)], [mean(combinedPsthSuccessMatrix(:,:,p)) + errSuccess, fliplr(mean(combinedPsthSuccessMatrix(:,:,p)) - errSuccess)], [0, 0, 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        title(['Mean PSTH (Successful Trials) with error bars - Fiber ', num2str(p), '- N = ', num2str(size(combinedRTToneSuccessful, 2)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
    end
    % Segment the RT data into thirds, fifths, tenths, etc.
        
    FastRTMissThreshold = 1.9; %Threshold for near misses (green), in seconds
    SlowRTMissThreshold = 0.1; %Threshold for false starts (magenta), in seconds
    sortedRT = sort(combinedRTToneSuccessful);

    % %UNCOMMENT BELOW IF YOU WANT TO COMPARE FASTEST VS SLOWEST THIRD, FIFTH, TENTH, ETC.
    % RTdivision = 3; %divide the data into thirds
    % idxFastest = floor(length(combinedRTToneSuccessful) / RTdivision);
    % idxSlowest = floor(length(combinedRTToneSuccessful) - idxFastest);
    % fastestEnd = sortedRT(idxFastest);
    % slowestStart = sortedRT(idxSlowest);
    % fastestStart = 0;
    % slowestEnd = 0.5; %this should be the RT threshold. Keep at 0.5
    
    %UNCOMMENT BELOW IF YOU WANT TO COMPARE FASTEST VS SLOWEST AT SPECIFIC CUTOFFS
    fastestEnd = 0.15;
    fastestStart = 0;
    slowestStart = fastestEnd;
    slowestEnd = 0.5;

    % Indices for the fastest and slowest 
    fastestIndices = combinedRTToneSuccessful >= fastestStart & combinedRTToneSuccessful < fastestEnd;
    slowestIndices = combinedRTToneSuccessful >= slowestStart & combinedRTToneSuccessful <= slowestEnd;

    % Extract PSTH data for the fastest and slowest 
    psthFastest = combinedPsthSuccessMatrix(fastestIndices, :, :);
    psthSlowest = combinedPsthSuccessMatrix(slowestIndices, :, :);
    
    % Segment the unsuccessful RT data based on threshold values
    
    % Indices for the fast and slow misses
    fastMissIndices = combinedRTToneUnsuccessful > FastRTMissThreshold;
    slowMissIndices = combinedRTToneUnsuccessful < SlowRTMissThreshold;
        
    % Extract PSTH data for the fast and slow misses
    psthFastestMiss = combinedPsthWithoutSecondMatrix(fastMissIndices, :, :);
    psthSlowestMiss = combinedPsthWithoutSecondMatrix(slowMissIndices, :, :);
    
    figure(5); 
    for p = 1:numfibers
        subplot(numfibers, 1, p);
        
        % Compute the mean and error for fastest third data
        meanFastest = mean(psthFastest(:,:,p), 1);
        errFastest = std(psthFastest(:,:,p), 0, 1) / sqrt(size(psthFastest, 1));
        
        % Compute the mean and error for slowest third data
        meanSlowest = mean(psthSlowest(:,:,p),1);
        errSlowest = std(psthSlowest(:,:,p),0,1) / sqrt(size(psthSlowest, 1));

        % Compute the mean and error for fastest missed trials
        meanFastestMiss = mean(psthFastestMiss(:,:,p),1);
        errFastestMiss = std(psthFastestMiss(:,:,p),0,1) / sqrt(size(psthFastestMiss, 1));
        
        % Compute the mean and error for slowest missed trials
        meanSlowestMiss = mean(psthSlowestMiss(:,:,p),1);
        errSlowestMiss = std(psthSlowestMiss(:,:,p),0,1) / sqrt(size(psthSlowestMiss, 1));

        % Plotting fastest third data
        plot(psthTime, meanFastest, 'b');
        hold on;
        if size(psthFastest, 1) > 1 % Only plot error bands if there's more than one trial
            fill([psthTime, fliplr(psthTime)], [meanFastest + errFastest, fliplr(meanFastest - errFastest)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        % Plotting slowest third data
        plot(psthTime, meanSlowest, 'r');
        if size(psthSlowest, 1) > 1 % Only plot error bands if there's more than one trial
            fill([psthTime, fliplr(psthTime)], [meanSlowest + errSlowest, fliplr(meanSlowest - errSlowest)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        xlabel('Time (s)');
        ylabel('Response');
        title(['Overlay PSTH - Fastest RTs (Blue, N = ', num2str(size(psthFastest, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(fastestIndices))), '), Slowest RTs (Red, N = ', num2str(size(psthSlowest, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(slowestIndices))), ') - Fiber ', num2str(p)]);
        plot([0 0], [-1 2.5], 'k-', 'LineWidth', 2);
        plot([2 2], [-1 2.5], 'k--', 'LineWidth', 2);
        plot([2.5 2.5], [-1 2.5], 'k--', 'LineWidth', 2);
        box off;
        ylim([-0.75 1.75]);  % Set y-axis limits here

        hold off;
    end


    figure(6); 
    for p = 1:numfibers
        subplot(numfibers, 1, p);
        
        % Compute the mean and error for fastest third data
        meanFastest = mean(psthFastest(:,:,p), 1);
        errFastest = std(psthFastest(:,:,p), 0, 1) / sqrt(size(psthFastest, 1));
        
        % Compute the mean and error for slowest third data
        meanSlowest = mean(psthSlowest(:,:,p),1);
        errSlowest = std(psthSlowest(:,:,p),0,1) / sqrt(size(psthSlowest, 1));

        % Compute the mean and error for fastest missed trials
        meanFastestMiss = mean(psthFastestMiss(:,:,p),1);
        errFastestMiss = std(psthFastestMiss(:,:,p),0,1) / sqrt(size(psthFastestMiss, 1));
        
        % Compute the mean and error for slowest missed trials
        meanSlowestMiss = mean(psthSlowestMiss(:,:,p),1);
        errSlowestMiss = std(psthSlowestMiss(:,:,p),0,1) / sqrt(size(psthSlowestMiss, 1));

        % Plotting fastest data
        plot(psthTime, meanFastest, 'b');
        hold on;
        if size(psthFastest, 1) > 1 % Only plot error bands if there's more than one trial
            fill([psthTime, fliplr(psthTime)], [meanFastest + errFastest, fliplr(meanFastest - errFastest)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        % Plotting slowest data
        plot(psthTime, meanSlowest, 'r');
        if size(psthSlowest, 1) > 1 % Only plot error bands if there's more than one trial
            fill([psthTime, fliplr(psthTime)], [meanSlowest + errSlowest, fliplr(meanSlowest - errSlowest)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        % Plotting fastest missed data
        plot(psthTime, meanFastestMiss, 'g');
        hold on;
        if size(psthFastestMiss, 1) > 1 
            fill([psthTime, fliplr(psthTime)], [meanFastestMiss + errFastestMiss, fliplr(meanFastestMiss - errFastestMiss)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    
        % Plotting slowest missed data
        plot(psthTime, meanSlowestMiss, 'm');
        if size(psthSlowestMiss, 1) > 1 
            fill([psthTime, fliplr(psthTime)], [meanSlowestMiss + errSlowestMiss, fliplr(meanSlowestMiss - errSlowestMiss)], 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        xlabel('Time (s)');
        ylabel('Response');
        title(['Overlay PSTH - Fast RTs (Blue, N = ', num2str(size(psthFastest, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(fastestIndices))), '), Slow RTs (Red, N = ', num2str(size(psthSlowest, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(slowestIndices))), '), Near Misses (Green, N = ', num2str(size(psthFastestMiss, 1)), ', Mean RT = ', num2str(mean(combinedRTToneUnsuccessful(fastMissIndices))), '), False Starts (Magenta, N = ', num2str(size(psthSlowestMiss, 1)), ', Mean RT = ', num2str(mean(combinedRTToneUnsuccessful(slowMissIndices))), ') - Fiber ', num2str(p)]);
        plot([0 0], [-1 2.5], 'k-', 'LineWidth', 2);
        plot([2 2], [-1 2.5], 'k--', 'LineWidth', 2);
        plot([2.5 2.5], [-1 2.5], 'k--', 'LineWidth', 2);
        ylim([-0.75 1.75]);  % Set y-axis limits here

        hold off;
        box off;
    end

    % Extracting max peaks for time window (-3.5, -2.5)
    [~, idxStart] = min(abs(psthTime - (-3.5)));
    [~, idxEnd] = min(abs(psthTime - (-2.5)));
    for f = 1:numfibers
        for t = 1:length(combinedRTToneSuccessful)
            currentData = combinedPsthSuccessMatrix(t, idxStart:idxEnd, f);
            maxPeaks(t, f, 1) = max(currentData);
        end
    end
    
    % Scatter plot for max peak vs RT
    figure(7);
    for f = 1:numfibers
        subplot(numfibers, 1, f);
        
        % Scatter plot
        scatter(combinedRTToneSuccessful, maxPeaks(:, f, 1));
        
        % Linear regression
        lm = fitlm(combinedRTToneSuccessful, maxPeaks(:, f, 1));
        
        % Sort RTToneSuccessful for proper plotting of regression and confidence intervals
        [sortedRT, sortedIndices] = sort(combinedRTToneSuccessful);
        
        % Get predicted values and confidence intervals based on sorted RT
        [ypred, yci] = predict(lm, sortedRT(:));
        
        % Plot regression line
        hold on;
        plot(sortedRT, ypred, 'r', 'LineWidth', 2);  % regression line
        
        % Plot confidence intervals
        plot(sortedRT, yci(:,1), 'r--');  % lower bound
        plot(sortedRT, yci(:,2), 'r--');  % upper bound
        
        hold off;
    
        xlabel('Reaction Time (RT)');
        ylabel('Max Peak');
        title(['Scatter plot of Max Peak vs RT (timeframe: -3.5s to -2.5s) - Fiber ', num2str(f), ' - N = ', num2str(size(maxPeaks, 1)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
        grid on;
    end
    
    % Extracting AUC for time window (-4, 0)
    [~, idxStart] = min(abs(psthTime - (-4)));
    [~, idxEnd] = min(abs(psthTime - (-0)));
    AUC = zeros(length(combinedRTToneSuccessful), numfibers);
    for f = 1:numfibers
        for t = 1:length(combinedRTToneSuccessful)
            currentData = combinedPsthSuccessMatrix(t, idxStart:idxEnd, f);
            AUC(t, f) = trapz(currentData);  % Compute the AUC
        end
    end
    
    % Scatter plot for AUC vs RT
    figure(8);
    for f = 1:numfibers
        subplot(numfibers, 1, f);
        
        % Scatter plot
        scatter(combinedRTToneSuccessful, AUC(:, f));
        
        % Linear regression
        lm = fitlm(combinedRTToneSuccessful, AUC(:, f));
        
        % Sort RTToneSuccessful for proper plotting of regression and confidence intervals
        [sortedRT, sortedIndices] = sort(combinedRTToneSuccessful);
        
        % Get predicted values and confidence intervals based on sorted RT
        [ypred, yci] = predict(lm, sortedRT(:));
        
        % Plot regression line
        hold on;
        plot(sortedRT, ypred, 'r', 'LineWidth', 2);  % regression line
        
        % Plot confidence intervals
        plot(sortedRT, yci(:,1), 'r--');  % lower bound
        plot(sortedRT, yci(:,2), 'r--');  % upper bound
        
        hold off;
    
        xlabel('Reaction Time (RT)');
        ylabel('AUC');
        title(['Scatter plot of AUC vs RT (timeframe: -3.5s to -2.5s)  - Fiber ', num2str(f), ' - N = ', num2str(size(AUC, 1)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
        grid on;
    end

    % Extract the experiment numbers from a sample mouse (assuming all mice have the same experiments)
    sampleMouse = mice{1};
    experimentNames = fieldnames(percentsuccessfultrials.(sampleMouse));
    
    % Convert experiment names to numbers for x-axis
    experimentNumbers = str2double(extractAfter(experimentNames, "experiment_"));
    
    % Create a figure for plotting
    figure(9);
        colors = lines(length(mice));
        hold on;
        for mouseIdx = 1:length(mice)
            mouseName = mice{mouseIdx};
            successPercentages = zeros(size(experimentNames));
            
            for expIdx = 1:length(experimentNames)
                successPercentages(expIdx) = percentsuccessfultrials.(mouseName).(experimentNames{expIdx});
            end
            
            plot(experimentNumbers, successPercentages, 'o-', 'Color', colors(mouseIdx, :), 'DisplayName', mouseName);
        end
        xticks(experimentNumbers);
        xlabel('Experiment Number');
        ylabel('Percentage of Successful Trials');
        legend('show', 'Location', 'best');
        title('Percentage of Successful Trials per Experiment for Different Mice');
        grid on;
        hold off;
        
        % Define rolling window parameters
        windowSize = 0.5; % e.g., 0.5 seconds
        stepSize = 0.1; % e.g., 0.1 seconds, for overlapping windows
        
        % Start and end times based on psthTime range
        startTime = min(psthTime);
        endTime = max(psthTime) - windowSize;
        
        % Initialize rolling window centers
        rollingWindowCenters = [];
        combinedDataTable = table(); % Collect data for all fibers
        
        % Initialize structures to save mixed-effects models for each fiber
        lmeAucModels = struct();
        lmeMaxPeakModels = struct();
        
        % Initialize arrays to store AUC and Max Peak difference values for each rolling window and fiber.
        numWindows = floor((endTime - startTime) / stepSize) + 1;
        rollingAucDifference = zeros(numfibers, numWindows);
        rollingMaxPeakDifference = zeros(numfibers, numWindows);
        rollingAucPValue = zeros(numfibers, numWindows);
        rollingMaxPeakPValue = zeros(numfibers, numWindows);
        
        windowIdx = 0; % Initialize window index
        
        % Slide the window from start to end time
        for windowStart = startTime:stepSize:endTime
            windowIdx = windowIdx + 1; % Increment window index
            windowEnd = windowStart + windowSize;
            rollingWindowCenters(end+1) = (windowStart + windowEnd) / 2; % Center of the current window
        
            % Find indices for the current window
            [~, idxStart] = min(abs(psthTime - windowStart));
            [~, idxEnd] = min(abs(psthTime - windowEnd));
        
            for a = 1:numfibers
                aucValuesAll = [];
                maxPeakValuesAll = [];
                fiberRTsAll = [];
                mouseAll = [];
                fiberNumAll = []; 
                idxStartAll = [];
                idxEndAll = [];
        
                for m = 1:length(mice)
                    mouseName = mice{m};
        
                    % Extract data for the current fiber and mouse
                    allDataForMouse = combinedPsthSuccess_struct.(mouseName);
                    fiberData = squeeze(allDataForMouse(:,:,a)); 
                    fiberRTs = combinedRTToneSuccessful_struct.(mouseName);
        
                    % Compute AUC and max peak for the current window for each trial
                    aucValues = trapz(fiberData(:, idxStart:idxEnd), 2);
                    maxPeakValues = max(fiberData(:, idxStart:idxEnd), [], 2);
        
                    aucValuesAll = [aucValuesAll; aucValues];
                    maxPeakValuesAll = [maxPeakValuesAll; maxPeakValues];
                    fiberRTsAll = [fiberRTsAll, fiberRTs];
                    mouseAll = [mouseAll; repmat(mouseName, length(aucValues), 1)];
                    fiberNumAll = [fiberNumAll; repmat(a, length(aucValues), 1)];
                    idxStartAll = [idxStartAll; repmat(idxStart, length(aucValues), 1)];
                    idxEndAll = [idxEndAll; repmat(idxEnd, length(aucValues), 1)];
                end
        
                fiberRTsAll = fiberRTsAll';
        
                % Populate the dataTable for the current fiber
                dataTable = table(idxStartAll, idxEndAll, aucValuesAll, maxPeakValuesAll, fiberRTsAll, mouseAll, fiberNumAll, 'VariableNames', {'StartIndex', 'EndIndex','AUC', 'MaxPeak', 'RT', 'Mouse', 'Fiber'});
                combinedDataTable = [combinedDataTable; dataTable]; % Concatenate the data
        
                % Create a mixed-effects model using the data
                % lmeAuc = fitlme(dataTable, 'AUC ~ 1 + RT + (1|Mouse)');
                % lmeMaxPeak = fitlme(dataTable, 'MaxPeak ~ 1 + RT + (1|Mouse)');
        
                % Create a simple linear regression model using the data
                lmeAuc = fitlm(dataTable, 'RT ~ AUC');
                lmeMaxPeak = fitlm(dataTable, 'RT ~ MaxPeak');
        
                % Save the models in the structures
                fiberName = strcat('Fiber', num2str(a));
                lmeAucModels.(fiberName) = lmeAuc;
                lmeMaxPeakModels.(fiberName) = lmeMaxPeak;
        
                % Extract the coefficients for the fixed effect of RT on AUC and MaxPeak
                coefAuc = lmeAuc.Coefficients.Estimate(2);
                coefMaxPeak = lmeMaxPeak.Coefficients.Estimate(2); 
                pValueAuc = lmeAuc.Coefficients.pValue(2);
                pValueMaxPeak = lmeMaxPeak.Coefficients.pValue(2);
        
                % Store them
                rollingAucDifference(a, windowIdx) = coefAuc;
                rollingMaxPeakDifference(a, windowIdx) = coefMaxPeak;
                rollingAucPValue(a, windowIdx) = pValueAuc;
                rollingMaxPeakPValue(a, windowIdx) = pValueMaxPeak;
            end
        end
        
        % Plotting the data for each fiber
        for a = 1:numfibers
            figure(a+10); % creates a new figure for each fiber
        
            % Subplot 1 for AUC Difference of the current fiber
            subplot(2, 1, 1);
            plot(rollingWindowCenters, rollingAucDifference(a, :), 'b');
            xlabel('Time (s)');
            ylabel('AUC Estimate');
            title(['Rolling AUC Estimate for Fiber ', num2str(a), ' - N = ', num2str(size(fiberRTsAll, 1)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
            grid on;
        
            % Subplot 2 for Max Peak Difference of the current fiber
            subplot(2, 1, 2);
            plot(rollingWindowCenters, rollingMaxPeakDifference(a, :), 'b');
            xlabel('Time (s)');
            ylabel('Max Peak Estimate');
            title(['Rolling Max Peak Estimate for Fiber ', num2str(a), ' - N = ', num2str(size(fiberRTsAll, 1)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
            grid on;
        end
        
        % Plotting the p-values for each fiber
        for a = 1:numfibers
            figure(a+12); % creates a new figure for each fiber
        
            % Subplot 1 for AUC p-value of the current fiber
            subplot(2, 1, 1);
            plot(rollingWindowCenters, rollingAucPValue(a, :), 'r'); % Using red color for clarity
            xlabel('Time (s)');
            ylabel('AUC p-value');
            title(['Rolling AUC p-value for Fiber ', num2str(a), ' - N = ', num2str(size(fiberRTsAll, 1)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
            grid on;
            ylim([0 0.2]);
        
            % Subplot 2 for Max Peak p-value of the current fiber
            subplot(2, 1, 2);
            plot(rollingWindowCenters, rollingMaxPeakPValue(a, :), 'r');
            xlabel('Time (s)');
            ylabel('Max Peak p-value');
            title(['Rolling Max Peak p-value for Fiber ', num2str(a), ' - N = ', num2str(size(fiberRTsAll, 1)), ' - Mean RT = ', num2str(mean(combinedRTToneSuccessful))]);
            grid on;
            ylim([0 0.2]);
        end

    
    % ITI-based segmentation
    over16IndicesSuccessful = combinedITISuccessful >= ITIlongthreshold;
    below16IndicesSuccessful = combinedITISuccessful <= ITIshortthreshold;
    fastestover16Indices = fastestIndices & over16IndicesSuccessful;
    fastestbelow16Indices = fastestIndices & below16IndicesSuccessful;
    slowestover16Indices = slowestIndices & over16IndicesSuccessful;
    slowestbelow16Indices = slowestIndices & below16IndicesSuccessful;
    
    over16IndicesUnsuccessful = combinedITIUnsuccessful >= ITIlongthreshold;
    below16IndicesUnsuccessful = combinedITIUnsuccessful <= ITIshortthreshold;
    fastMissover16Indices = fastMissIndices & over16IndicesUnsuccessful;
    fastMissbelow16Indices = fastMissIndices & below16IndicesUnsuccessful;
    slowMissover16Indices = slowMissIndices & over16IndicesUnsuccessful;
    slowMissbelow16Indices = slowMissIndices & below16IndicesUnsuccessful;

    % Extracting PSTH based on ITI
    psthFastestOver = combinedPsthSuccessMatrix(fastestover16Indices, :, :);
    psthFastestBelow = combinedPsthSuccessMatrix(fastestbelow16Indices, :, :);
    psthSlowestOver = combinedPsthSuccessMatrix(slowestover16Indices, :, :);
    psthSlowestBelow = combinedPsthSuccessMatrix(slowestbelow16Indices, :, :);

    psthFastestMissOver = combinedPsthWithoutSecondMatrix(fastMissover16Indices, :, :);
    psthFastestMissBelow = combinedPsthWithoutSecondMatrix(fastMissbelow16Indices, :, :);
    psthSlowestMissOver = combinedPsthWithoutSecondMatrix(slowMissover16Indices, :, :);
    psthSlowestMissBelow = combinedPsthWithoutSecondMatrix(slowMissbelow16Indices, :, :);

    %Calculate Mean/SD ITI during fast/slow RTs
    mean_ITI_fastRT = mean(combinedITISuccessful( :, fastestIndices));
    sd_ITI_fastRT = std(combinedITISuccessful( :, fastestIndices));
    mean_ITI_slowRT = mean(combinedITISuccessful( :, slowestIndices));
    sd_ITI_slowRT = std(combinedITISuccessful( :, slowestIndices));
    ITIttest = struct();
    [ITIttest.h, ITIttest.p] = ttest2(combinedITISuccessful( :, fastestIndices), combinedITISuccessful( :, slowestIndices), 'VarType', 'unequal');

    figure(15);
    for p = 1:numfibers
        % For ITI >= ITIlongthreshold
        subplot(numfibers, 2, 2*p-1);

        meanFastest = mean(psthFastestOver(:,:,p), 1);
        errFastest = std(psthFastestOver(:,:,p), 0, 1) / sqrt(size(psthFastestOver, 1));
        plot(psthTime, meanFastest, 'b');
        hold on;
        fill([psthTime, fliplr(psthTime)], [meanFastest + errFastest, fliplr(meanFastest - errFastest)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        meanSlowest = mean(psthSlowestOver(:,:,p), 1);
        errSlowest = std(psthSlowestOver(:,:,p), 0, 1) / sqrt(size(psthSlowestOver, 1));
        plot(psthTime, meanSlowest, 'r');
        fill([psthTime, fliplr(psthTime)], [meanSlowest + errSlowest, fliplr(meanSlowest - errSlowest)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        meanFastestMiss = mean(psthFastestMissOver(:,:,p), 1);
        errFastestMiss = std(psthFastestMissOver(:,:,p), 0, 1) / sqrt(size(psthFastestMissOver, 1));
        plot(psthTime, meanFastestMiss, 'g');
        fill([psthTime, fliplr(psthTime)], [meanFastestMiss + errFastestMiss, fliplr(meanFastestMiss - errFastestMiss)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        meanSlowestMiss = mean(psthSlowestMissOver(:,:,p), 1);
        errSlowestMiss = std(psthSlowestMissOver(:,:,p), 0, 1) / sqrt(size(psthSlowestMissOver, 1));
        plot(psthTime, meanSlowestMiss, 'm');
        fill([psthTime, fliplr(psthTime)], [meanSlowestMiss + errSlowestMiss, fliplr(meanSlowestMiss - errSlowestMiss)], 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        line11 = ['Fiber ', num2str(p), ' - ITI >=', num2str(ITIlongthreshold), 's - Fastest RTs (Blue, N = ', num2str(size(psthFastestOver, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(fastestover16Indices))), '), Slowest RTs (Red, N = ', num2str(size(psthSlowestOver, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(slowestover16Indices))), ')'];
        line12 = ['Fast Misses (Green, N = ', num2str(size(psthFastestMissOver, 1)), ', Mean RT = ', num2str(mean(combinedRTToneUnsuccessful(fastMissover16Indices))), '), Slow Misses (Magenta, N = ', num2str(size(psthSlowestMissOver, 1)), ', Mean RT = ', num2str(mean(combinedRTToneUnsuccessful(slowMissover16Indices))), ')'];

        title({line11, line12});

        % For ITI <= ITIshortthreshold
        subplot(numfibers, 2, 2*p);

        meanFastest = mean(psthFastestBelow(:,:,p), 1);
        errFastest = std(psthFastestBelow(:,:,p), 0, 1) / sqrt(size(psthFastestBelow, 1));
        plot(psthTime, meanFastest, 'b');
        hold on;
        fill([psthTime, fliplr(psthTime)], [meanFastest + errFastest, fliplr(meanFastest - errFastest)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        meanSlowest = mean(psthSlowestBelow(:,:,p), 1);
        errSlowest = std(psthSlowestBelow(:,:,p), 0, 1) / sqrt(size(psthSlowestBelow, 1));
        plot(psthTime, meanSlowest, 'r');
        fill([psthTime, fliplr(psthTime)], [meanSlowest + errSlowest, fliplr(meanSlowest - errSlowest)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        meanFastestMiss = mean(psthFastestMissBelow(:,:,p), 1);
        errFastestMiss = std(psthFastestMissBelow(:,:,p), 0, 1) / sqrt(size(psthFastestMissBelow, 1));
        plot(psthTime, meanFastestMiss, 'g');
        fill([psthTime, fliplr(psthTime)], [meanFastestMiss + errFastestMiss, fliplr(meanFastestMiss - errFastestMiss)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        meanSlowestMiss = mean(psthSlowestMissBelow(:,:,p), 1);
        errSlowestMiss = std(psthSlowestMissBelow(:,:,p), 0, 1) / sqrt(size(psthSlowestMissBelow, 1));
        plot(psthTime, meanSlowestMiss, 'm');
        fill([psthTime, fliplr(psthTime)], [meanSlowestMiss + errSlowestMiss, fliplr(meanSlowestMiss - errSlowestMiss)], 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

        line21 = ['Fiber ', num2str(p), ' - ITI <=', num2str(ITIshortthreshold), 's - Fastest RTs (Blue, N = ', num2str(size(psthFastestBelow, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(fastestbelow16Indices))), '), Slowest RTs (Red, N = ', num2str(size(psthSlowestBelow, 1)), ', Mean RT = ', num2str(mean(combinedRTToneSuccessful(slowestbelow16Indices))), ')'];
        line22 = ['Fast Misses (Green, N = ', num2str(size(psthFastestMissBelow, 1)), ', Mean RT = ', num2str(mean(combinedRTToneUnsuccessful(fastMissbelow16Indices))), '), Slow Misses (Magenta, N = ', num2str(size(psthSlowestMissBelow, 1)), ', Mean RT = ', num2str(mean(combinedRTToneUnsuccessful(slowMissbelow16Indices))), ')'];

        title({line21, line22});

    end    
    
    % Adjust the RTs as specified
    combinedRTSuccessful_ms = combinedRTSuccessful * 1000;
    combinedRTUnsuccessful_ms = -1 * ((combinedRTUnsuccessful * 1000) - 2000);
    combinedRTs_ms = [combinedRTSuccessful_ms, -combinedRTUnsuccessful_ms];
    combinedRTs_ms_sorted = sort(combinedRTs_ms);
    
    % plot the combined RTs
    figure (16);
        histogram(combinedRTs_ms, 'BinWidth', 100);
        xlim([-5000, 5000]);
        title('Histogram of Reaction Times');
        xlabel('Reaction Time (milliseconds)');
        ylabel('Frequency');
        grid on;
        
        hold on;
        line([0 0], ylim, 'Color', 'red', 'LineWidth', 2);
        text(0, max(ylim) * 0.9, 'Second Tone', 'HorizontalAlignment', 'center', 'Color', 'black', 'FontWeight', 'bold');
        line([-2000 -2000], ylim, 'Color', 'red', 'LineWidth', 2);
        text(-2000, max(ylim) * 0.9, 'First Tone', 'HorizontalAlignment', 'center', 'Color', 'black', 'FontWeight', 'bold');
        hold off;

    % Find the mode of the unsuccessful reaction times
    mode_RTToneUnsuccessful = mode(combinedRTToneUnsuccessful);
    
    % Define the mode window of ±0.1 s
    mode_window = 0.1;
    
    % Filter the RTs to find those within the mode window
    RTs_in_mode_window = combinedRTToneUnsuccessful(combinedRTToneUnsuccessful >= (mode_RTToneUnsuccessful - mode_window) & combinedRTToneUnsuccessful <= (combinedRTToneUnsuccessful + mode_window));
    
    % Sort the RTs within the mode window
    sorted_RTs_in_mode_window = sort(RTs_in_mode_window);
    
    % Determine the indices for the first and last thirds within this filtered, sorted array
    numModeTrials = length(sorted_RTs_in_mode_window);
    idxFirstThird = floor(numModeTrials / 3);
    idxLastThird = numModeTrials - floor(numModeTrials / 3);
    
    % Get the cutoff reaction time values for the first and last thirds
    firstThirdCutoff = sorted_RTs_in_mode_window(idxFirstThird);
    firstThirdStart = mode_RTToneUnsuccessful - mode_window;
    lastThirdCutoff = sorted_RTs_in_mode_window(idxLastThird);
    lastThirdStart = mode_RTToneUnsuccessful + mode_window;
    
    % Find the indices in the original unsorted RT array that correspond to the first and last third cutoffs
    firstThirdIndices = combinedRTToneUnsuccessful > firstThirdStart & combinedRTToneUnsuccessful < firstThirdCutoff;
    lastThirdIndices = (combinedRTToneUnsuccessful > lastThirdStart) & (combinedRTToneUnsuccessful < lastThirdCutoff);
    
    % Extract the PSTH data for these trials
    psthFirstThird = combinedPsthWithoutSecondMatrix(firstThirdIndices, :, :);
    psthLastThird = combinedPsthWithoutSecondMatrix(lastThirdIndices, :, :);

    figure(17); 
    for p = 1:numfibers
        subplot(numfibers, 1, p);

        % Compute the mean and error for the first third data
        meanFirstThird = mean(psthFirstThird(:,:,p), 1);
        errFirstThird = std(psthFirstThird(:,:,p), 0, 1) / sqrt(size(psthFirstThird, 1));

        % Compute the mean and error for the last third data
        meanLastThird = mean(psthLastThird(:,:,p),1);
        errLastThird = std(psthLastThird(:,:,p),0,1) / sqrt(size(psthLastThird, 1));

        % Plotting first third data
        plot(psthTime, meanFirstThird, 'b');
        hold on;
        if size(psthFirstThird, 1) > 1 % Only plot error bands if there's more than one trial
            fill([psthTime, fliplr(psthTime)], [meanFirstThird + errFirstThird, fliplr(meanFirstThird - errFirstThird)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        % Plotting last third data
        plot(psthTime, meanLastThird, 'r');
        if size(psthLastThird, 1) > 1 % Only plot error bands if there's more than one trial
            fill([psthTime, fliplr(psthTime)], [meanLastThird + errLastThird, fliplr(meanLastThird - errLastThird)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        xlabel('Time (s)');
        ylabel('Response');
        title(['Overlay PSTH - First Third RTs +/- 100 ms from mode (Blue, N = ', num2str(size(psthFirstThird, 1)), ...
               '), Last Third RTs +/- 100 ms from mode (Red, N = ', num2str(size(psthLastThird, 1)), ...
               ') - Fiber ', num2str(p)]);
        plot([0 0], ylim, 'k-', 'LineWidth', 2); % Adjust the ylim based on your data range
        plot([2 2], ylim, 'k--', 'LineWidth', 2);
        plot([2.5 2.5], ylim, 'k--', 'LineWidth', 2);
        xlim([min(psthTime) max(psthTime)]); % Set the x-axis limits to match the psthTime range

        hold off;
    end

    % Calculate the indices for each segment
    sortedRT = sort(combinedRTToneSuccessful,'ascend'); % Sort the reaction times in ascending order
    numTrials = length(sortedRT);
    segments = [10, 9, 8, 7, 6, 5, 4, 3,2]; % Define the segments: half, third, quarter, etc.
    segmentIndices = arrayfun(@(x) floor(numTrials/x), segments);

    % Calculate mean RT for each segment
    meanRTs = arrayfun(@(x) mean(sortedRT(1:x)), segmentIndices);
    
    % Define the color gradient
    numSegments = length(segments);
    colors = jet(numSegments);
    
    figure(18); % Create a new figure for the segments visualization
    for p = 1:numfibers
        subplot(numfibers, 1, p); hold on;
        for s = 1:numSegments
            % Calculate the index range for the current segment
            idxEnd = segmentIndices(s);
            currentSegmentIndices = sortedRT(1:idxEnd);
            
            % Find the trials that fall within the current segment
            segmentRTindices = ismember(combinedRTToneSuccessful, currentSegmentIndices);
    
            % Compute the mean and error for the current segment
            psthCurrentSegment = combinedPsthSuccessMatrix(segmentRTindices, :, p);
            meanCurrentSegment = mean(psthCurrentSegment, 1);
            errCurrentSegment = std(psthCurrentSegment, 0, 1) / sqrt(size(psthCurrentSegment, 1));
            
            % Plotting current segment data with gradient color
            plot(psthTime, meanCurrentSegment, 'color', colors(s,:), 'LineWidth', 1.5);
            
            % Creating the shaded error area
            fill([psthTime, fliplr(psthTime)], ...
                 [meanCurrentSegment + errCurrentSegment, fliplr(meanCurrentSegment - errCurrentSegment)], ...
                 colors(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        xlabel('Time (s)');
        ylabel('Response');
        title(['Gradient PSTH by Fastest Segments - Fiber ', num2str(p)]);
        plot([0 0], [-1 2.5], 'k-', 'LineWidth', 2);
        plot([2 2], [-1 2.5], 'k--', 'LineWidth', 2);
        plot([2.5 2.5], [-1 2.5], 'k--', 'LineWidth', 2);
        ylim([-0.75 1.75]);
        hold off;
    end

    % Define minimum number of observations per bin
    minObservations = 100;
    
    % Initialize variables
    binStart = 0;
    binEnd = 0.02; % start with 20 ms bins
    binEdges = [0]; % to store the actual used bin edges
    binCounts = []; % to store the number of observations per bin
    
    % Determine the bins based on the minimum number of observations
    while binEnd <= 0.5
        RTindices = combinedRTToneSuccessful >= binStart & combinedRTToneSuccessful < binEnd;
        binCount = sum(RTindices);
        
        if binCount < minObservations
            % Extend the bin until we have enough observations
            binEnd = binEnd + 0.02; % Extend the bin by 20ms
        else
            % We have enough observations, so we finalize this bin
            binEdges = [binEdges, binEnd];
            binCounts = [binCounts, binCount];
            binStart = binEnd; % Move to the next bin
            binEnd = binStart + 0.02; % Reset the bin end to start 20ms after the new start
        end
    end
    
    % If the last bin has fewer than the minimum, merge it with the previous bin
    if binCounts(end) < minObservations && length(binCounts) > 1
        binEdges(end-1) = binEdges(end);
        binEdges(end) = [];
        binCounts(end-1) = binCounts(end-1) + binCounts(end);
        binCounts(end) = [];
    end
    
    % determine the color gradient based on the number of bins we have
    numBins = length(binEdges) - 1;
    colors = jet(numBins);
    
    % Create the plots as before, using the new binEdges
    figure(19);
    for p = 1:numfibers
        subplot(numfibers, 1, p); hold on;
        for binIdx = 1:numBins
            binStart = binEdges(binIdx);
            binEnd = binEdges(binIdx+1);
            RTindices = combinedRTToneSuccessful >= binStart & combinedRTToneSuccessful < binEnd;
            
            psthCurrentBin = combinedPsthSuccessMatrix(RTindices, :, p);
            meanCurrentBin = mean(psthCurrentBin, 1);
            errCurrentBin = std(psthCurrentBin, 0, 1) / sqrt(size(psthCurrentBin, 1));
    
            plot(psthTime, meanCurrentBin, 'color', colors(binIdx,:), 'LineWidth', 1.5);
            fill([psthTime, fliplr(psthTime)], ...
                 [meanCurrentBin + errCurrentBin, fliplr(meanCurrentBin - errCurrentBin)], ...
                 colors(binIdx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        xlabel('Time (s)');
        ylabel('Response');
        title(['Gradient PSTH by RT bins - Fiber ', num2str(p)]);
        plot([0 0], [-1 2.5], 'k-', 'LineWidth', 2);
        plot([2 2], [-1 2.5], 'k--', 'LineWidth', 2);
        plot([2.5 2.5], [-1 2.5], 'k--', 'LineWidth', 2);
        ylim([-0.75 1.75]);
        hold off;
    end


    % Calculate the indices for each segment for the slowest trials
    sortedRT = sort(combinedRTToneSuccessful, 'descend'); % Sort the reaction times in descending order
    numTrials = length(sortedRT);
    segments = [10, 9, 8, 7, 6, 5, 4, 3, 2]; % Define the segments: half, third, quarter, etc.
    segmentIndices = arrayfun(@(x) floor(numTrials/x), segments);
    
    % Calculate mean RT for each segment
    meanRTs = arrayfun(@(x) mean(sortedRT(1:x)), segmentIndices);
    
    % Define the color gradient
    numSegments = length(segments);
    colors = flip(jet(numSegments));
    
    figure(20); % Create a new figure for the segments visualization
    for p = 1:numfibers
        subplot(numfibers, 1, p); hold on;
        for s = 1:numSegments
            % Calculate the index range for the current segment
            idxEnd = segmentIndices(s);
            currentSegmentIndices = sortedRT(1:idxEnd);
            
            % Find the trials that fall within the current segment
            segmentRTindices = ismember(combinedRTToneSuccessful, currentSegmentIndices);
    
            % Compute the mean and error for the current segment
            psthCurrentSegment = combinedPsthSuccessMatrix(segmentRTindices, :, p);
            meanCurrentSegment = mean(psthCurrentSegment, 1);
            errCurrentSegment = std(psthCurrentSegment, 0, 1) / sqrt(size(psthCurrentSegment, 1));
            
            % Plotting current segment data with gradient color
            plot(psthTime, meanCurrentSegment, 'color', colors(s,:), 'LineWidth', 1.5);
            
            % Creating the shaded error area
            fill([psthTime, fliplr(psthTime)], ...
                 [meanCurrentSegment + errCurrentSegment, fliplr(meanCurrentSegment - errCurrentSegment)], ...
                 colors(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
        
        xlabel('Time (s)');
        ylabel('Response');
        title(['Gradient PSTH by Slowest Segments - Fiber ', num2str(p)]);
        plot([0 0], [-1 2.5], 'k-', 'LineWidth', 2);
        plot([2 2], [-1 2.5], 'k--', 'LineWidth', 2);
        plot([2.5 2.5], [-1 2.5], 'k--', 'LineWidth', 2);
        ylim([-1 2.5]);
        hold off;
    end


    % Time windows of interest
    timeWindows = {[-4,-2],[-2, 0],[0, 2],[2,2.5]};
    
    % Compute various metrics for multiple time windows
    numWindows = length(timeWindows);
    aucValues = NaN(length(combinedRTToneSuccessful), numfibers, numWindows);
    maxPeaks = NaN(length(combinedRTToneSuccessful), numfibers, numWindows);
    maxSlopes = NaN(length(combinedRTToneSuccessful), numfibers, numWindows);
    timeOfMaxPeak = NaN(length(combinedRTToneSuccessful), numfibers, numWindows); % Initialize with NaNs

    for w = 1:numWindows
        currentWindow = timeWindows{w};
        
        % Find indices for the start and end of the current time window
        [~, idxStart] = min(abs(psthTime - currentWindow(1)));
        [~, idxEnd] = min(abs(psthTime - currentWindow(2)));
        
        for a = 1:numfibers
            dataWindow = combinedPsthSuccessMatrix(:, idxStart:idxEnd, a);
            
            % AUC
            aucValues(:, a, w) = trapz(dataWindow, 2);
            
            % Max Peak
            maxPeaks(:, a, w) = max(dataWindow, [], 2);
            
            % Max Slope
            slopes = diff(dataWindow, 1, 2) / (psthTime(2) - psthTime(1));
            maxSlopes(:, a, w) = max(slopes, [], 2);
            
            % Time of Max Peak relative to first tone
            [~, idxMaxPeak] = max(dataWindow, [], 2);
            timeOfMaxPeak(:, a, w) = psthTime(idxStart + idxMaxPeak - 1);
        end
    end
    
    % Correlation computation for each metric
    correlationResultsAUC = NaN(numfibers, numWindows, 2);
    correlationResultsMax = NaN(numfibers, numWindows, 2);
    correlationResultsSlope = NaN(numfibers, numWindows, 2);
    correlationResultsTimeToPeak = NaN(numfibers, numWindows, 2); 
    
    for w = 1:numWindows
        for a = 1:numfibers
            [R_AUC, P_AUC] = corrcoef(aucValues(:, a, w), combinedRTToneSuccessful');
            correlationResultsAUC(a, w, 1) = R_AUC(1, 2);  % AUC correlation coefficient
            correlationResultsAUC(a, w, 2) = P_AUC(1, 2);  % AUC p-value
            
            [R_Max, P_Max] = corrcoef(maxPeaks(:, a, w), combinedRTToneSuccessful');
            correlationResultsMax(a, w, 1) = R_Max(1, 2);  % Max Peak correlation coefficient
            correlationResultsMax(a, w, 2) = P_Max(1, 2);  % Max Peak p-value
            
            [R_Slope, P_Slope] = corrcoef(maxSlopes(:, a, w), combinedRTToneSuccessful');
            correlationResultsSlope(a, w, 1) = R_Slope(1, 2);  % Max Slope correlation coefficient
            correlationResultsSlope(a, w, 2) = P_Slope(1, 2);  % Max Slope p-value
            
            if timeWindows{w}(1) < 0 && timeWindows{w}(2) <= 0
                [R_TimeToPeak, P_TimeToPeak] = corrcoef(timeOfMaxPeak(:, a, w), combinedRTToneSuccessful');
                correlationResultsTimeToPeak(a, w, 1) = R_TimeToPeak(1, 2);
                correlationResultsTimeToPeak(a, w, 2) = P_TimeToPeak(1, 2);
            end
        end
    end
    
    % Display results for each fiber and time window
    for a = 1:numfibers
        fprintf('Fiber %d:\n', a);
        for w = 1:numWindows
            fprintf('Time Window: %d to %d seconds\n', timeWindows{w}(1), timeWindows{w}(2));
            fprintf('Correlation between AUC and RT: %.3f (p = %.3f)\n', correlationResultsAUC(a, w, 1), correlationResultsAUC(a, w, 2));
            fprintf('Correlation between Max Peak and RT: %.3f (p = %.3f)\n', correlationResultsMax(a, w, 1), correlationResultsMax(a, w, 2));
            fprintf('Correlation between Max Slope and RT: %.3f (p = %.3f)\n', correlationResultsSlope(a, w, 1), correlationResultsSlope(a, w, 2));
            
            if timeWindows{w}(1) < 0
                fprintf('Correlation between Time of Max Peak and RT: %.3f (p = %.3f)\n', correlationResultsTimeToPeak(a, w, 1), correlationResultsTimeToPeak(a, w, 2));
            else
                fprintf('Time of Max Peak and its correlation with RT is not applicable for this window.\n');
            end
            fprintf('\n');
    
            % Display average Max Peak and Time of Max Peak for the window
            avgMaxPeak = mean(maxPeaks(:, a, w));
            sdMaxPeak = std(maxPeaks(:, a, w));
            fprintf('Average Max Peak in this window: %.3f ± %.3f\n', avgMaxPeak, sdMaxPeak);
            
            avgTimeOfMaxPeak = mean(timeOfMaxPeak(:, a, w));
            sdTimeOfMaxPeak = std(timeOfMaxPeak(:, a, w));
            fprintf('Average Time of Max Peak in this window: %.3f ± %.3f seconds\n', avgTimeOfMaxPeak, sdTimeOfMaxPeak);
            fprintf('\n');
        end
    end

    % Compute various metrics for multiple time windows
    numWindows = length(timeWindows);
    
    numfibers = size(psthFastest, 3);
    tValues = NaN(numfibers, numWindows);
    pValues = NaN(numfibers, numWindows);
    
    for w = 1:numWindows
        currentWindow = timeWindows{w};
        
        % Find indices for the start and end of the current time window
        [~, idxStart] = min(abs(psthTime - currentWindow(1)));
        [~, idxEnd] = min(abs(psthTime - currentWindow(2)));
        
        for a = 1:numfibers
            dataFast = squeeze(psthFastest(:, idxStart:idxEnd, a));
            dataSlow = squeeze(psthSlowest(:, idxStart:idxEnd, a));
    
            % Flatten the data to a single vector for each condition
            dataFast = dataFast(:);
            dataSlow = dataSlow(:);
    
            % Perform Welch's t-test for unequal variance between dataFast and dataSlow
            RTttest = struct();
            [RTttest.h, RTttest.p] = ttest2(dataFast, dataSlow, 'VarType', 'unequal');
            
        end
    end
    
    %stats and graphs for manuscript figure 6/10/24
    %Panel c scatter plot
    % Define minimum number of observations per bin
    minObservations = 100;
    
    % Initialize variables
    binStart = 0;
    binEnd = 0.02; % start with 20 ms bins
    binEdges = [0]; % to store the actual used bin edges
    binCounts = []; % to store the number of observations per bin
    
    % Determine the bins based on the minimum number of observations
    while binEnd <= 0.5
        RTindices = combinedRTToneSuccessful >= binStart & combinedRTToneSuccessful < binEnd;
        binCount = sum(RTindices);
    
        if binCount < minObservations
            % Extend the bin until we have enough observations
            binEnd = binEnd + 0.02; % Extend the bin by 20ms
        else
            % We have enough observations, so we finalize this bin
            binEdges = [binEdges, binEnd];
            binCounts = [binCounts, binCount];
            binStart = binEnd; % Move to the next bin
            binEnd = binStart + 0.02; % Reset the bin end to start 20ms after the new start
        end
    end
    
    % If the last bin has fewer than the minimum, merge it with the previous bin
    if binCounts(end) < minObservations && length(binCounts) > 1
        binEdges(end-1) = binEdges(end);
        binEdges(end) = [];
        binCounts(end-1) = binCounts(end-1) + binCounts(end);
        binCounts(end) = [];
    end
    
    % determine the color gradient based on the number of bins we have
    numBins = length(binEdges) - 1;
    colors = jet(numBins);
    
    % Initialize arrays to store RT bins and corresponding mean z-scores from 1 to 2s
    RTBins = zeros(1, numBins);
    meanZScores = zeros(1, numBins);
    sems = zeros(1, numBins);
    
    % Loop through bins to calculate the mean z-score from 1 to 2s for each bin
    for binIdx = 1:numBins
        binStart = binEdges(binIdx);
        binEnd = binEdges(binIdx+1);
        RTindices = combinedRTToneSuccessful >= binStart & combinedRTToneSuccessful < binEnd;
    
        psthCurrentBin = combinedPsthSuccessMatrix(RTindices, :, :);
        timeIndices = (psthTime >= 0) & (psthTime <= 2);
        dataInRange = psthCurrentBin(:, timeIndices, :);
    
        meanZScore = mean(dataInRange(:));
        semZScore = std(dataInRange(:)) / sqrt(length(dataInRange(:)));
    
        RTBins(binIdx) = (binStart + binEnd) / 2; % Middle of the bin
        meanZScores(binIdx) = meanZScore;
        sems(binIdx) = semZScore;
    end
    
    % Create scatter plot
    figure;
    for binIdx = 1:numBins
        errorbar(RTBins(binIdx), meanZScores(binIdx), sems(binIdx), 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(binIdx,:), 'Color', colors(binIdx,:));
        hold on;
    end
    
    % Perform linear regression
    [p, S] = polyfit(RTBins, meanZScores, 1);
    [yfit, delta] = polyval(p, RTBins, S);
    plot(RTBins, yfit, '-r');
    
    % Calculate the confidence interval
    yresid = meanZScores - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(meanZScores)-1) * var(meanZScores);
    rsq = 1 - SSresid/SStotal;
    n = length(meanZScores);
    alpha = 0.05;
    pred_ci = tinv(1-alpha/2, n-2) * sqrt(SSresid/(n-2) * (1/n + (RTBins - mean(RTBins)).^2 / sum((RTBins - mean(RTBins)).^2)));
    
    % Plot the confidence interval
    plot(RTBins, yfit + pred_ci, 'b--');
    plot(RTBins, yfit - pred_ci, 'b--');
    
    % Add R^2 and p-value to the plot
    text(0.3, max(meanZScores) + 0.05, sprintf('R^2 = %.2f\np = %.3f', rsq, pval), 'FontSize', 12, 'HorizontalAlignment', 'center');
    
    xlabel('Mean RT (s)');
    ylabel('Mean z-score from 1 to 2s');
    title('Scatter plot of RT bins vs. mean z-score from 1 to 2s');
    
    % Adjust the axis limits if necessary
    xlim([0, max(RTBins) + 0.01]);
    ylim([min(meanZScores) - 0.05, max(meanZScores) + 0.05]);
    
    grid off;
    hold off;
    box off;

    %Panel d scatter plot
    % Calculate the indices for each segment
    sortedRT = sort(combinedRTToneSuccessful, 'ascend'); % Sort the reaction times in ascending order
    numTrials = length(sortedRT);
    segments = [10, 9, 8, 7, 6, 5, 4, 3, 2]; % Define the segments: half, third, quarter, etc.
    segmentIndices = arrayfun(@(x) floor(numTrials/x), segments);
    
    % Calculate mean RT for each segment
    meanRTs = arrayfun(@(x) mean(sortedRT(1:x)), segmentIndices);
    
    % Define the color gradient
    numSegments = length(segments);
    colors = jet(numSegments);
    
    % Initialize arrays to store the means and SEMs for each segment
    means = zeros(1, length(segmentIndices));
    sems = zeros(1, length(segmentIndices));
    
    % Time range for mean calculation
    timeRange = -4:0;
    
    % Loop through each segment to calculate the mean and SEM for the specified time range
    for s = 1:length(segmentIndices)
        idxEnd = segmentIndices(s);
        currentSegmentIndices = sortedRT(1:idxEnd);
    
        % Find the trials that fall within the current segment
        segmentRTindices = ismember(combinedRTToneSuccessful, currentSegmentIndices);
    
        % Compute the mean and SEM for the specified time range for the current segment
        psthCurrentSegment = combinedPsthSuccessMatrix(segmentRTindices, :, :);
        timeIndices = ismember(psthTime, timeRange);
        dataInRange = psthCurrentSegment(:, timeIndices, :);
    
        meanValue = mean(dataInRange(:));
        semValue = std(dataInRange(:)) / sqrt(length(dataInRange(:)));
    
        means(s) = meanValue;
        sems(s) = semValue;
    end
    
    % Create scatter plot
    figure;
    for s = 1:length(segmentIndices)
        errorbar(meanRTs(s), means(s), sems(s), 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(s,:), 'Color', colors(s,:));
        hold on;
    end
    
    % Perform quadratic regression
    [p, S] = polyfit(meanRTs, means, 2);
    [yfit, delta] = polyval(p, meanRTs, S);
    plot(meanRTs, yfit, '-r');
    
    % Calculate the confidence interval
    yresid = means - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(means)-1) * var(means);
    rsq = 1 - SSresid/SStotal;
    n = length(means);
    alpha = 0.05;
    pred_ci = tinv(1-alpha/2, n-2) * sqrt(SSresid/(n-2) * (1/n + (meanRTs - mean(meanRTs)).^2 / sum((meanRTs - mean(meanRTs)).^2)));
    
    % Plot the confidence interval
    plot(meanRTs, yfit + pred_ci, 'b--');
    plot(meanRTs, yfit - pred_ci, 'b--');
    
    % Calculate p-value for R²
    [pval, ~] = anova1(means, meanRTs, 'off');
    
    % Add statistics to the plot
    if pval < 0.001
        pval_text = 'p < 0.001';
    else
        pval_text = sprintf('p = %.3f', pval);
    end
    
    text(mean(meanRTs)+ 0.04, max(means), sprintf('R^2 = %.2f\n%s', rsq, pval_text), 'FontSize', 12);
    
    xlabel('Mean RT (s)');
    ylabel('Mean z-score from -4 to 0s');
    title('Scatter plot of RT bins vs. mean z-score from -4 to 0s');
    
    % Adjust the axis limits if necessary
    xlim([0, 0.1]);
    ylim([- 0.5, -0.25]);
    
    grid off;
    hold off;
    box off;



    %panel e barplot 1

    % Define RT ranges for panel E
    fastestStart = 0;
    fastestEnd = 0.02;
    slowestStart = 0.02;
    slowestEnd = 0.5;
    
    % Indices for the fastest and slowest 
    fastestIndicesE = combinedRTToneSuccessful >= fastestStart & combinedRTToneSuccessful < fastestEnd;
    slowestIndicesE = combinedRTToneSuccessful >= slowestStart & combinedRTToneSuccessful <= slowestEnd;
    
    % Extract z-score data for the fastest and slowest
    zScoresFastestE = combinedPsthSuccessMatrix(fastestIndicesE, :, :);
    zScoresSlowestE = combinedPsthSuccessMatrix(slowestIndicesE, :, :);
    
    % Define the time window for analysis (-4 to 0 seconds)
    timeWindowIndicesE = (psthTime >= -4) & (psthTime < 0);
    
    % Extract z-scores for the defined time window
    zScoresFastestWindowE = zScoresFastestE(:, timeWindowIndicesE, 1);
    zScoresSlowestWindowE = zScoresSlowestE(:, timeWindowIndicesE, 1);
    
    % Compute the mean z-score for the time window for each trial
    meanZScoresFastestE = mean(zScoresFastestWindowE, 2);
    meanZScoresSlowestE = mean(zScoresSlowestWindowE, 2);
    
    % Normalize the z-scores as difference from -0.5
    normZScoresFastestE = meanZScoresFastestE - (-0.5);
    normZScoresSlowestE = meanZScoresSlowestE - (-0.5);
    
    % Compute mean and SEM for the normalized z-scores
    meanNormZScoreFastestE = mean(normZScoresFastestE);
    meanNormZScoreSlowestE = mean(normZScoresSlowestE);
    semNormZScoreFastestE = std(normZScoresFastestE) / sqrt(length(normZScoresFastestE));
    semNormZScoreSlowestE = std(normZScoresSlowestE) / sqrt(length(normZScoresSlowestE));
    
    % Perform a two-sample t-test
    [~, pvalE] = ttest2(normZScoresFastestE, normZScoresSlowestE);
    
    % Create a bar plot
    figure;
    barDataE = [meanNormZScoreFastestE, meanNormZScoreSlowestE];
    semDataE = [semNormZScoreFastestE, semNormZScoreSlowestE];
    barHandleE = bar(barDataE);
    hold on;
    errorbar(1:2, barDataE, semDataE, 'k', 'linestyle', 'none');
    set(gca, 'XTickLabel', {'RT 0-20ms', 'RT 20-500ms'});
    ylabel('Normalized Mean Z-Score');
    title('Normalized Mean Z-Score Over -4 to 0s (Panel E)');
    if pvalE < 0.001
        pval_textE = 'p < 0.001';
    else
        pval_textE = sprintf('p = %.3f', pvalE);
    end
    text(1.5, max(barDataE) + 0.05, pval_textE, 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Set colors
    barHandleE.FaceColor = 'flat';
    barHandleE.CData(1, :) = [0 0 1]; % Blue for fastest RTs
    barHandleE.CData(2, :) = [1 0 0]; % Red for slowest RTs
    
    % Customize plot appearance
    box off;
    grid off;
    hold off;

    % Define RT ranges for panel E barplot 2
    fastestStart = 0;
    fastestEnd = 0.02;
    slowestStart = 0.02;
    slowestEnd = 0.5;
    
    % Indices for the fastest and slowest
    fastestIndicesE = combinedRTToneSuccessful >= fastestStart & combinedRTToneSuccessful < fastestEnd;
    slowestIndicesE = combinedRTToneSuccessful >= slowestStart & combinedRTToneSuccessful <= slowestEnd;
    
    % Extract z-score data for the fastest and slowest
    zScoresFastestE = combinedPsthSuccessMatrix(fastestIndicesE, :, :);
    zScoresSlowestE = combinedPsthSuccessMatrix(slowestIndicesE, :, :);
    
    % Define the time window for analysis (0 to 2 seconds)
    timeWindowIndicesE = (psthTime >= 0) & (psthTime < 2);
    
    % Extract z-scores for the defined time window
    zScoresFastestWindowE = zScoresFastestE(:, timeWindowIndicesE, 1);
    zScoresSlowestWindowE = zScoresSlowestE(:, timeWindowIndicesE, 1);
    
    % Normalize the z-scores as difference from -0.5
    normZScoresFastestE = zScoresFastestWindowE - (-0.5);
    normZScoresSlowestE = zScoresSlowestWindowE - (-0.5);
    
    % Compute mean and SEM for the normalized z-scores
    meanNormZScoreFastestE = mean(normZScoresFastestE(:));  % Flatten to 1D, then take mean
    meanNormZScoreSlowestE = mean(normZScoresSlowestE(:));  % Flatten to 1D, then take mean
    
    semNormZScoreFastestE = std(normZScoresFastestE(:)) / sqrt(length(normZScoresFastestE(:)));
    semNormZScoreSlowestE = std(normZScoresSlowestE(:)) / sqrt(length(normZScoresSlowestE(:)));
    
    % Perform a two-sample t-test
    [~, pvalE] = ttest2(normZScoresFastestE(:), normZScoresSlowestE(:));
    
    % Create a bar plot
    figure;
    barDataE = [meanNormZScoreFastestE, meanNormZScoreSlowestE];
    semDataE = [semNormZScoreFastestE, semNormZScoreSlowestE];
    
    barHandleE = bar(barDataE);  % Create the bar chart with the means
    hold on;
    
    % Add error bars
    errorbar(1:2, barDataE, semDataE, 'k', 'linestyle', 'none');
    
    % Customize axis labels and title
    set(gca, 'XTickLabel', {'RT 0-20ms', 'RT 20-500ms'});
    ylabel('Normalized Mean Z-Score');
    title('Normalized Mean Z-Score Over 0s to 2s (Panel E)');
    
    % Add p-value text
    if pvalE < 0.001
        pval_textE = 'p < 0.001';
    else
        pval_textE = sprintf('p = %.3f', pvalE);
    end
    
    % Display p-value text on the plot
    text(1.5, max(barDataE) + 0.05, pval_textE, 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Set colors for the bars
    barHandleE.FaceColor = 'flat';
    barHandleE.CData(1, :) = [0 0 1]; % Blue for fastest RTs
    barHandleE.CData(2, :) = [1 0 0]; % Red for slowest RTs
    
    % Customize plot appearance
    box off;
    grid off;
    hold off;

    %Panel e barplot 3: Define RT ranges 
    % --> AKA PANEL J IN UPDATED FIGUERE 1

    fastestStart_e = 0;
    fastestEnd_e = 0.02;
    slowestStart_e = fastestEnd_e;
    slowestEnd_e = 0.5;
    
    % Panel e: Indices for the fastest and slowest trials
    fastestIndices_e = combinedRTToneSuccessful >= fastestStart_e & combinedRTToneSuccessful < fastestEnd_e;
    slowestIndices_e = combinedRTToneSuccessful > slowestStart_e & combinedRTToneSuccessful <= slowestEnd_e;
    
    % Panel e: Extract z-score data for the fastest and slowest trials
    zScoresFastest_e = combinedPsthSuccessMatrix(fastestIndices_e, :, :);
    zScoresSlowest_e = combinedPsthSuccessMatrix(slowestIndices_e, :, :);
    
    % Panel e: Define the time windows for analysis
    timeWindow1Indices_e = (psthTime >= -4) & (psthTime < 0);
    timeWindow2Indices_e = (psthTime >= 0) & (psthTime <= 2);
    
    % Panel e: Extract z-scores for the defined time windows
    zScoresFastestWindow1_e = zScoresFastest_e(:, timeWindow1Indices_e, 1);
    zScoresSlowestWindow1_e = zScoresSlowest_e(:, timeWindow1Indices_e, 1);
    
    zScoresFastestWindow2_e = zScoresFastest_e(:, timeWindow2Indices_e, 1);
    zScoresSlowestWindow2_e = zScoresSlowest_e(:, timeWindow2Indices_e, 1);
    
    % Panel e: Compute the mean z-scores for the time windows
    meanZScoresFastestWindow1_e = mean(zScoresFastestWindow1_e(:));
    meanZScoresSlowestWindow1_e = mean(zScoresSlowestWindow1_e(:));
    
    meanZScoresFastestWindow2_e = mean(zScoresFastestWindow2_e(:));
    meanZScoresSlowestWindow2_e = mean(zScoresSlowestWindow2_e(:));
    
    % Panel e: Calculate the differences in z-scores between fast and slow trials for each time window
    diffWindow1_e = meanZScoresFastestWindow1_e - meanZScoresSlowestWindow1_e;
    diffWindow2_e = meanZScoresFastestWindow2_e - meanZScoresSlowestWindow2_e;
    
    % Normalize the z-scores as difference from -0.5
    normZScoresFastestE1 = zScoresFastestWindow1_e - (0.5);
    normZScoresSlowestE1 = zScoresSlowestWindow1_e - (0.5);
    
    normZScoresFastestE2 = zScoresFastestWindow2_e - (0.5);
    normZScoresSlowestE2 = zScoresSlowestWindow2_e - (0.5);
    
    % Perform a two-sample t-test
    [~, pvalE1] = ttest2(normZScoresFastestE1(:), normZScoresSlowestE1(:));
    [~, pvalE2] = ttest2(normZScoresFastestE2(:), normZScoresSlowestE2(:));
    
    % Panel e: Create the bar plot
    figure;
    barData = [diffWindow1_e, diffWindow2_e];
    bar([1, 2], barData, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    
    % Panel e: Add error bars (SEM)
    semFastestWindow1_e = std(zScoresFastestWindow1_e(:)) / sqrt(length(zScoresFastestWindow1_e(:)));
    semSlowestWindow1_e = std(zScoresSlowestWindow1_e(:)) / sqrt(length(zScoresSlowestWindow1_e(:)));
    semDiffWindow1_e = sqrt(semFastestWindow1_e^2 + semSlowestWindow1_e^2);
    
    semFastestWindow2_e = std(zScoresFastestWindow2_e(:)) / sqrt(length(zScoresFastestWindow2_e(:)));
    semSlowestWindow2_e = std(zScoresSlowestWindow2_e(:)) / sqrt(length(zScoresSlowestWindow2_e(:)));
    semDiffWindow2_e = sqrt(semFastestWindow2_e^2 + semSlowestWindow2_e^2);
    
    errorbar([1, 2], barData, [semDiffWindow1_e, semDiffWindow2_e], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Set the y-axis limits
    ylim([-0.1 0.2]); % Limit the y-axis range from -0.1 to 0.2
    
    % Panel e: Customize plot appearance
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'-4 to 0 s', '0 to 2 s'});
    ylabel('Difference in Mean Z-Score');
    title('Panel e: Difference in Z-Scores Between Fast and Slow Trials');
    box off;
    grid off;
    
    % Adjust the height where the p-values are displayed (reduce the offset)
    pval_offset = 0.02;  % Set a smaller offset for better visibility
    
    % Add p-values to the plot above the bars
    if pvalE1 < 0.001
        pvalTextE1 = 'p < 0.001';
    else
        pvalTextE1 = sprintf('p = %.3f', pvalE1);
    end
    
    if pvalE2 < 0.001
        pvalTextE2 = 'p < 0.001';
    else
        pvalTextE2 = sprintf('p = %.3f', pvalE2);
    end
    
    % Display p-values above the corresponding bars with reduced offset
    text(1, barData(1) + semDiffWindow1_e + pval_offset, pvalTextE1, 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(2, barData(2) + semDiffWindow2_e + pval_offset, pvalTextE2, 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Hold off to complete the plot
    hold off;
    
    % Panel e: Display the actual values for verification
    disp(['Panel e - Difference in Z-Score for -4 to 0 s: ', num2str(diffWindow1_e)]);
    disp(['Panel e - Difference in Z-Score for 0 to 2 s: ', num2str(diffWindow2_e)]);
    disp(['Panel e - p-value for -4 to 0 s: ', num2str(pvalE1)]);
    disp(['Panel e - p-value for 0 to 2 s: ', num2str(pvalE2)]);

    %panel f barplot 1
    % Define RT ranges for panel F
    RTdivision = 3; % divide the data into thirds
    idxFastestF = floor(length(combinedRTToneSuccessful) / RTdivision);
    idxSlowestF = floor(length(combinedRTToneSuccessful) - idxFastestF);
    sortedRT = sort(combinedRTToneSuccessful);
    fastestEndF = sortedRT(idxFastestF);
    slowestStartF = sortedRT(idxSlowestF);
    fastestStartF = 0;
    slowestEndF = 0.5;
    
    % Indices for the fastest and slowest 
    fastestIndicesF = combinedRTToneSuccessful >= fastestStartF & combinedRTToneSuccessful < fastestEndF;
    slowestIndicesF = combinedRTToneSuccessful >= slowestStartF & combinedRTToneSuccessful <= slowestEndF;
    
    % Extract z-score data for the fastest and slowest
    zScoresFastestF = combinedPsthSuccessMatrix(fastestIndicesF, :, :);
    zScoresSlowestF = combinedPsthSuccessMatrix(slowestIndicesF, :, :);
    
    % Define the time window for analysis (-4 to 0 seconds)
    timeWindowIndicesF = (psthTime >= -4) & (psthTime < 0);
    
    % Extract z-scores for the defined time window
    zScoresFastestWindowF = zScoresFastestF(:, timeWindowIndicesF, 1);
    zScoresSlowestWindowF = zScoresSlowestF(:, timeWindowIndicesF, 1);
    
    % Compute the mean z-score for the time window for each trial
    meanZScoresFastestF = mean(zScoresFastestWindowF, 2);
    meanZScoresSlowestF = mean(zScoresSlowestWindowF, 2);
    
    % Normalize the z-scores as difference from -0.5
    normZScoresFastestF = meanZScoresFastestF - (-0.5);
    normZScoresSlowestF = meanZScoresSlowestF - (-0.5);
    
    % Compute mean and SEM for the normalized z-scores
    meanNormZScoreFastestF = mean(normZScoresFastestF);
    meanNormZScoreSlowestF = mean(normZScoresSlowestF);
    semNormZScoreFastestF = std(normZScoresFastestF) / sqrt(length(normZScoresFastestF));
    semNormZScoreSlowestF = std(normZScoresSlowestF) / sqrt(length(normZScoresSlowestF));
    
    % Perform a two-sample t-test
    [~, pvalF] = ttest2(normZScoresFastestF, normZScoresSlowestF);
    
    % Create a bar plot
    figure;
    barDataF = [meanNormZScoreFastestF, meanNormZScoreSlowestF];
    semDataF = [semNormZScoreFastestF, semNormZScoreSlowestF];
    barHandleF = bar(barDataF);
    hold on;
    errorbar(1:2, barDataF, semDataF, 'k', 'linestyle', 'none');
    set(gca, 'XTickLabel', {'Fastest RTs', 'Slowest RTs'});
    ylabel('Normalized Mean Z-Score');
    title('Normalized Mean Z-Score Over -4 to 0s (Panel F)');
    if pvalF < 0.001
        pval_textF = 'p < 0.001';
    else
        pval_textF = sprintf('p = %.3f', pvalF);
    end
    text(1.5, max(barDataF) + 0.05, pval_textF, 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Set colors
    barHandleF.FaceColor = 'flat';
    barHandleF.CData(1, :) = [0 0 1]; % Blue for fastest RTs
    barHandleF.CData(2, :) = [0 1 1]; % Cyan for slowest RTs
    
    % Customize plot appearance
    box off;
    grid off;
    hold off;

    %panel f barplot 2
    % Define RT ranges for panel F
    RTdivision = 3; % divide the data into thirds
    idxFastestF = floor(length(combinedRTToneSuccessful) / RTdivision);
    idxSlowestF = floor(length(combinedRTToneSuccessful) - idxFastestF);
    sortedRT = sort(combinedRTToneSuccessful);
    fastestEndF = sortedRT(idxFastestF);
    slowestStartF = sortedRT(idxSlowestF);
    fastestStartF = 0;
    slowestEndF = 0.5;
    
    % Indices for the fastest and slowest 
    fastestIndicesF = combinedRTToneSuccessful >= fastestStartF & combinedRTToneSuccessful < fastestEndF;
    slowestIndicesF = combinedRTToneSuccessful >= slowestStartF & combinedRTToneSuccessful <= slowestEndF;
    
    % Extract z-score data for the fastest and slowest
    zScoresFastestF = combinedPsthSuccessMatrix(fastestIndicesF, :, :);
    zScoresSlowestF = combinedPsthSuccessMatrix(slowestIndicesF, :, :);
   
    % Define the time window for analysis (0 to 2 seconds)
    timeWindowIndicesF = (psthTime >= 0) & (psthTime < 2);
    
    % Extract z-scores for the defined time window
    zScoresFastestWindowF = zScoresFastestF(:, timeWindowIndicesF, 1);
    zScoresSlowestWindowF = zScoresSlowestF(:, timeWindowIndicesF, 1);
    
    % Compute the mean z-score for the time window for each trial
    meanZScoresFastestF = mean(zScoresFastestWindowF, 2);
    meanZScoresSlowestF = mean(zScoresSlowestWindowF, 2);
    
    % Normalize the z-scores as difference from -0.5
    normZScoresFastestF = meanZScoresFastestF - (-0.5);
    normZScoresSlowestF = meanZScoresSlowestF - (-0.5);
    
    % Compute mean and SEM for the normalized z-scores
    meanNormZScoreFastestF = mean(normZScoresFastestF);
    meanNormZScoreSlowestF = mean(normZScoresSlowestF);
    semNormZScoreFastestF = std(normZScoresFastestF) / sqrt(length(normZScoresFastestF));
    semNormZScoreSlowestF = std(normZScoresSlowestF) / sqrt(length(normZScoresSlowestF));
    
    % Perform a two-sample t-test
    [~, pvalF] = ttest2(normZScoresFastestF, normZScoresSlowestF);
    
    % Create a bar plot
    figure;
    barDataF = [meanNormZScoreFastestF, meanNormZScoreSlowestF];
    semDataF = [semNormZScoreFastestF, semNormZScoreSlowestF];
    barHandleF = bar(barDataF);
    hold on;
    errorbar(1:2, barDataF, semDataF, 'k', 'linestyle', 'none');
    set(gca, 'XTickLabel', {'Fastest RTs', 'Slowest RTs'});
    ylabel('Normalized Mean Z-Score');
    title('Normalized Mean Z-Score Over 0s to 2s (Panel F)');
    if pvalF < 0.001
        pval_textF = 'p < 0.001';
    else
        pval_textF = sprintf('p = %.3f', pvalF);
    end
    text(1.5, max(barDataF) + 0.05, pval_textF, 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Set colors
    barHandleF.FaceColor = 'flat';
    barHandleF.CData(1, :) = [0 0 1]; % Blue for fastest RTs
    barHandleF.CData(2, :) = [0 1 1]; % Cyan for slowest RTs
    
    % Customize plot appearance
    box off;
    grid off;
    hold off;

    %panel f barplot 3
    % --> AKA PANEL L iN UPDATED FIGURE 1

    % Panel f: Define RT ranges for fastest and slowest thirds
    sortedRT_f = sort(combinedRTToneSuccessful);
    RTdivision_f = 3; % Divide the data into thirds
    idxFastest_f = floor(length(combinedRTToneSuccessful) / RTdivision_f);
    idxSlowest_f = floor(length(combinedRTToneSuccessful) - idxFastest_f);
    fastestEnd_f = sortedRT_f(idxFastest_f);
    slowestStart_f = sortedRT_f(idxSlowest_f);
    fastestStart_f = 0;
    slowestEnd_f = 0.1; % Adjusted to <100ms
    
    % Panel f: Indices for the fastest and slowest trials
    fastestIndices_f = combinedRTToneSuccessful >= fastestStart_f & combinedRTToneSuccessful < fastestEnd_f;
    slowestIndices_f = combinedRTToneSuccessful >= slowestStart_f & combinedRTToneSuccessful <= slowestEnd_f;
    
    % Panel f: Extract z-score data for the fastest and slowest trials
    zScoresFastest_f = combinedPsthSuccessMatrix(fastestIndices_f, :, :);
    zScoresSlowest_f = combinedPsthSuccessMatrix(slowestIndices_f, :, :);
    
    % Panel f: Define the time windows for analysis
    timeWindow1Indices_f = (psthTime >= -4) & (psthTime < 0);
    timeWindow2Indices_f = (psthTime >= 0) & (psthTime <= 2);
    
    % Panel f: Extract z-scores for the defined time windows
    zScoresFastestWindow1_f = zScoresFastest_f(:, timeWindow1Indices_f, 1);
    zScoresSlowestWindow1_f = zScoresSlowest_f(:, timeWindow1Indices_f, 1);
    
    zScoresFastestWindow2_f = zScoresFastest_f(:, timeWindow2Indices_f, 1);
    zScoresSlowestWindow2_f = zScoresSlowest_f(:, timeWindow2Indices_f, 1);
    
    % Panel f: Compute the mean z-scores for the time windows
    meanZScoresFastestWindow1_f = mean(zScoresFastestWindow1_f(:));
    meanZScoresSlowestWindow1_f = mean(zScoresSlowestWindow1_f(:));
    
    meanZScoresFastestWindow2_f = mean(zScoresFastestWindow2_f(:));
    meanZScoresSlowestWindow2_f = mean(zScoresSlowestWindow2_f(:));
    
    % Panel f: Calculate the differences in z-scores between fast and slow trials for each time window
    diffWindow1_f = meanZScoresFastestWindow1_f - meanZScoresSlowestWindow1_f;
    diffWindow2_f = meanZScoresFastestWindow2_f - meanZScoresSlowestWindow2_f;
    
    % Normalize the z-scores as difference from -0.5
    normZScoresFastestF1 = zScoresFastestWindow1_f - (0.5);
    normZScoresSlowestF1 = zScoresSlowestWindow1_f - (0.5);
    
    normZScoresFastestF2 = zScoresFastestWindow2_f - (0.5);
    normZScoresSlowestF2 = zScoresSlowestWindow2_f - (0.5);
    
    % Perform a two-sample t-test
    [~, pvalF1] = ttest2(normZScoresFastestF1(:), normZScoresSlowestF1(:));
    [~, pvalF2] = ttest2(normZScoresFastestF2(:), normZScoresSlowestF2(:));
    
    % Panel f: Create the bar plot
    figure;
    barData_f = [diffWindow1_f, diffWindow2_f];
    bar([1, 2], barData_f, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    
    % Panel f: Add error bars (SEM)
    semFastestWindow1_f = std(zScoresFastestWindow1_f(:)) / sqrt(length(zScoresFastestWindow1_f(:)));
    semSlowestWindow1_f = std(zScoresSlowestWindow1_f(:)) / sqrt(length(zScoresSlowestWindow1_f(:)));
    semDiffWindow1_f = sqrt(semFastestWindow1_f^2 + semSlowestWindow1_f^2);
    
    semFastestWindow2_f = std(zScoresFastestWindow2_f(:)) / sqrt(length(zScoresFastestWindow2_f(:)));
    semSlowestWindow2_f = std(zScoresSlowestWindow2_f(:)) / sqrt(length(zScoresSlowestWindow2_f(:)));
    semDiffWindow2_f = sqrt(semFastestWindow2_f^2 + semSlowestWindow2_f^2);
    
    errorbar([1, 2], barData_f, [semDiffWindow1_f, semDiffWindow2_f], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Panel f: Customize plot appearance
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'-4 to 0 s', '0 to 2 s'});
    ylabel('Difference in Mean Z-Score');
    title('Panel f: Difference in Z-Scores Between Fast and Slow Trials');
    box off;
    grid off;
    
    % Add p-values to the plot above the bars
    pval_offset_f = 0.02; % Set a smaller offset for better visibility
    
    % Panel f: Add p-value text
    if pvalF1 < 0.001
        pvalTextF1 = 'p < 0.001';
    else
        pvalTextF1 = sprintf('p = %.3f', pvalF1);
    end
    
    if pvalF2 < 0.001
        pvalTextF2 = 'p < 0.001';
    else
        pvalTextF2 = sprintf('p = %.3f', pvalF2);
    end
    
    % Display p-values above the corresponding bars
    text(1, barData_f(1) + semDiffWindow1_f + pval_offset_f, pvalTextF1, 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(2, barData_f(2) + semDiffWindow2_f + pval_offset_f, pvalTextF2, 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % Hold off to complete the plot
    hold off;
    
    % Panel f: Display the actual values for verification
    disp(['Panel f - Difference in Z-Score for -4 to 0 s: ', num2str(diffWindow1_f)]);
    disp(['Panel f - Difference in Z-Score for 0 to 2 s: ', num2str(diffWindow2_f)]);
    disp(['Panel f - p-value for -4 to 0 s: ', num2str(pvalF1)]);
    disp(['Panel f - p-value for 0 to 2 s: ', num2str(pvalF2)]);

    %sample Plot
    % Define the indices for the second fastest and the second slowest trials
    sortedRTs = sort(combinedRTToneSuccessful);
    
    % Identify the second fastest and second slowest trials
    secondFastestRT = sortedRTs(2); % Second fastest reaction time
    secondSlowestRT = sortedRTs(end-1); % Second slowest reaction time
    
    secondFastestTrialIndex = find(combinedRTToneSuccessful == secondFastestRT, 1);
    secondSlowestTrialIndex = find(combinedRTToneSuccessful == secondSlowestRT, 1);
    
    % Extract the photometry data for the second fastest and second slowest trials
    secondFastestTrialData = combinedPsthSuccessMatrix(secondFastestTrialIndex, :, 1);
    secondSlowestTrialData = combinedPsthSuccessMatrix(secondSlowestTrialIndex, :, 1);
    
    % Create the figure
    figure;
    
    % Plot the photometry data for the second fastest and second slowest trials
    plot(psthTime, secondFastestTrialData, 'b', 'LineWidth', 1.5);
    hold on;
    plot(psthTime, secondSlowestTrialData, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Z-Score');
    title('Photometry Data for Second Fastest (Blue) and Second Slowest (Red) RT Trials');
    xlim([-10, 10]); % Adjust as necessary based on your data
    ylim([-2, 2.5]); % Adjust as necessary based on your data
    
    % Add vertical lines to indicate key time points
    plot([0 0], ylim, 'k-', 'LineWidth', 1);
    plot([2 2], ylim, 'k--', 'LineWidth', 1);
    plot([2.5 2.5], ylim, 'k--', 'LineWidth', 1);
    
    % Customize plot appearance
    legend({'Second Fastest RT', 'Second Slowest RT'}, 'Location', 'Best');
    box off;
    grid off;
    hold off;

    % Display the actual reaction times
    disp(['Second Fastest RT: ', num2str(secondFastestRT)]);
    disp(['Second Slowest RT: ', num2str(secondSlowestRT)]);

    %distribution of successful RTs
    % Plot Histogram of RTs
    figure;
    histogram(combinedRTToneSuccessful, 'BinWidth', 0.01);
    xlabel('Reaction Time (s)');
    ylabel('Count');
    title('Distribution of Successful RTs');
    box off;


    %panel b barplots
    % Define RT ranges
    fastestStart = 0;
    fastestEnd = 0.02;
    slowestStart = 0.02;
    slowestEnd = 0.5;
    
    % Indices for the fastest, slowest, near miss, and false start
    fastestIndices = combinedRTToneSuccessful >= fastestStart & combinedRTToneSuccessful < fastestEnd;
    slowestIndices = combinedRTToneSuccessful > slowestStart & combinedRTToneSuccessful <= slowestEnd;
    nearMissIndices = combinedRTToneUnsuccessful > FastRTMissThreshold;
    falseStartIndices = combinedRTToneUnsuccessful < SlowRTMissThreshold;
    
    % Extract z-score data for the fastest, slowest, near miss, and false start
    zScoresFastest_m = combinedPsthSuccessMatrix(fastestIndices, :, :);
    zScoresSlowest_m = combinedPsthSuccessMatrix(slowestIndices, :, :);
    zScoresFastestMiss_m = combinedPsthWithoutSecondMatrix(nearMissIndices, :, :);
    zScoresSlowestMiss_m = combinedPsthWithoutSecondMatrix(falseStartIndices, :, :);
    
    % Define the time windows for analysis
    timeWindow1Indices_m = (psthTime >= -4) & (psthTime < 0);
    timeWindow2Indices_m = (psthTime >= 0) & (psthTime <= 2);
    
    % Extract z-scores for the defined time windows
    zScoresFastestWindow1_m = zScoresFastest_m(:, timeWindow1Indices_m, 1);
    zScoresSlowestWindow1_m = zScoresSlowest_m(:, timeWindow1Indices_m, 1);
    zScoresFastestMissWindow1_m = zScoresFastestMiss_m(:, timeWindow1Indices_m, 1);
    zScoresSlowestMissWindow1_m = zScoresSlowestMiss_m(:, timeWindow1Indices_m, 1);
    
    zScoresFastestWindow2_m = zScoresFastest_m(:, timeWindow2Indices_m, 1);
    zScoresSlowestWindow2_m = zScoresSlowest_m(:, timeWindow2Indices_m, 1);
    zScoresFastestMissWindow2_m = zScoresFastestMiss_m(:, timeWindow2Indices_m, 1);
    zScoresSlowestMissWindow2_m = zScoresSlowestMiss_m(:, timeWindow2Indices_m, 1);
    
    % Normalize the z-scores for the time window (-4 to 0) by adding 1
    zScoresFastestWindow1_m = zScoresFastestWindow1_m + 1;
    zScoresSlowestWindow1_m = zScoresSlowestWindow1_m + 1;
    zScoresFastestMissWindow1_m = zScoresFastestMissWindow1_m + 1;
    zScoresSlowestMissWindow1_m = zScoresSlowestMissWindow1_m + 1;
    
    % Compute the mean z-scores for the time windows
    meanZScoresFastestWindow1_m = mean(zScoresFastestWindow1_m(:));
    meanZScoresSlowestWindow1_m = mean(zScoresSlowestWindow1_m(:));
    meanZScoresFastestMissWindow1_m = mean(zScoresFastestMissWindow1_m(:));
    meanZScoresSlowestMissWindow1_m = mean(zScoresSlowestMissWindow1_m(:));
    
    meanZScoresFastestWindow2_m = mean(zScoresFastestWindow2_m(:));
    meanZScoresSlowestWindow2_m = mean(zScoresSlowestWindow2_m(:));
    meanZScoresFastestMissWindow2_m = mean(zScoresFastestMissWindow2_m(:));
    meanZScoresSlowestMissWindow2_m = mean(zScoresSlowestMissWindow2_m(:));
    
    % Calculate the SEM for the time windows
    semFastestWindow1_m = std(zScoresFastestWindow1_m(:)) / sqrt(length(zScoresFastestWindow1_m(:)));
    semSlowestWindow1_m = std(zScoresSlowestWindow1_m(:)) / sqrt(length(zScoresSlowestWindow1_m(:)));
    semFastestMissWindow1_m = std(zScoresFastestMissWindow1_m(:)) / sqrt(length(zScoresFastestMissWindow1_m(:)));
    semSlowestMissWindow1_m = std(zScoresSlowestMissWindow1_m(:)) / sqrt(length(zScoresSlowestMissWindow1_m(:)));
    
    semFastestWindow2_m = std(zScoresFastestWindow2_m(:)) / sqrt(length(zScoresFastestWindow2_m(:)));
    semSlowestWindow2_m = std(zScoresSlowestWindow2_m(:)) / sqrt(length(zScoresSlowestWindow2_m(:)));
    semFastestMissWindow2_m = std(zScoresFastestWindow2_m(:)) / sqrt(length(zScoresFastestWindow2_m(:)));
    semSlowestMissWindow2_m = std(zScoresSlowestWindow2_m(:)) / sqrt(length(zScoresSlowestWindow2_m(:)));
    
    % Perform pairwise t-tests for time window (-4 to 0)
    [~, pval1_2_window1] = ttest2(zScoresFastestWindow1_m(:), zScoresSlowestWindow1_m(:));
    [~, pval1_3_window1] = ttest2(zScoresFastestWindow1_m(:), zScoresFastestMissWindow1_m(:));
    [~, pval1_4_window1] = ttest2(zScoresFastestWindow1_m(:), zScoresSlowestMissWindow1_m(:));
    [~, pval2_3_window1] = ttest2(zScoresSlowestWindow1_m(:), zScoresFastestMissWindow1_m(:));
    [~, pval2_4_window1] = ttest2(zScoresSlowestWindow1_m(:), zScoresSlowestMissWindow1_m(:));
    [~, pval3_4_window1] = ttest2(zScoresFastestMissWindow1_m(:), zScoresSlowestMissWindow1_m(:));
    
    % Perform pairwise t-tests for time window (0 to 2)
    [~, pval1_2_window2] = ttest2(zScoresFastestWindow2_m(:), zScoresSlowestWindow2_m(:));
    [~, pval1_3_window2] = ttest2(zScoresFastestWindow2_m(:), zScoresFastestMissWindow2_m(:));
    [~, pval1_4_window2] = ttest2(zScoresFastestWindow2_m(:), zScoresSlowestMissWindow2_m(:));
    [~, pval2_3_window2] = ttest2(zScoresSlowestWindow2_m(:), zScoresFastestMissWindow2_m(:));
    [~, pval2_4_window2] = ttest2(zScoresSlowestWindow2_m(:), zScoresSlowestMissWindow2_m(:));
    [~, pval3_4_window2] = ttest2(zScoresFastestMissWindow2_m(:), zScoresSlowestMissWindow2_m(:));
    
    % Create the bar plot for time window (-4 to 0)
    figure;
    subplot(1, 2, 1);
    bar(1:4, [meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresFastestMissWindow1_m, meanZScoresSlowestMissWindow1_m], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    errorbar(1:4, [meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresFastestMissWindow1_m, meanZScoresSlowestMissWindow1_m], ...
             [semFastestWindow1_m, semSlowestWindow1_m, semFastestMissWindow1_m, semSlowestMissWindow1_m], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Customize bar colors
    b = bar(1:4, [meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresFastestMissWindow1_m, meanZScoresSlowestMissWindow1_m], 'FaceColor', 'flat');
    b.CData(1,:) = [0 0 1];  % Blue
    b.CData(2,:) = [1 0 0];  % Red
    b.CData(3,:) = [0 1 0];  % Green
    b.CData(4,:) = [1 0 1];  % Magenta
    
    set(gca, 'XTickLabel', {'Fast RTs', 'Slow RTs', 'Near Misses', 'False Starts'});
    ylabel('Normalized Mean Z-Score');
    title('Mean Z-Score (-4 to 0 s)');
    box off;
    grid on;
    hold off;
    
    % Add p-values to the plot for time window (-4 to 0)
    text(1.5, median([meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresFastestMissWindow1_m, meanZScoresSlowestMissWindow1_m]) + 0.1, sprintf('p_{1,2} = %.3f\np_{1,3} = %.3f\np_{1,4} = %.3f\np_{2,3} = %.3f\np_{2,4} = %.3f\np_{3,4} = %.3f', pval1_2_window1, pval1_3_window1, pval1_4_window1, pval2_3_window1, pval2_4_window1, pval3_4_window1), 'FontSize', 12);
    
    % Create the bar plot for time window (0 to 2)
    subplot(1, 2, 2);
    bar(1:4, [meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresFastestMissWindow2_m, meanZScoresSlowestMissWindow2_m], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    errorbar(1:4, [meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresFastestMissWindow2_m, meanZScoresSlowestMissWindow2_m], ...
             [semFastestWindow2_m, semSlowestWindow2_m, semFastestMissWindow2_m, semSlowestMissWindow2_m], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Customize bar colors
    b = bar(1:4, [meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresFastestMissWindow2_m, meanZScoresSlowestMissWindow2_m], 'FaceColor', 'flat');
    b.CData(1,:) = [0 0 1];  % Blue
    b.CData(2,:) = [1 0 0];  % Red
    b.CData(3,:) = [0 1 0];  % Green
    b.CData(4,:) = [1 0 1];  % Magenta
    
    set(gca, 'XTickLabel', {'Fast RTs', 'Slow RTs', 'Near Misses', 'False Starts'});
    ylabel('Mean Z-Score');
    title('Mean Z-Score (0 to 2 s)');
    box off;
    grid on;
    hold off;
    
    % Add p-values to the plot for time window (0 to 2)
    text(1.5, median([meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresFastestMissWindow2_m, meanZScoresSlowestMissWindow2_m]) + 0.1, sprintf('p_{1,2} = %.3f\np_{1,3} = %.3f\np_{1,4} = %.3f\np_{2,3} = %.3f\np_{2,4} = %.3f\np_{3,4} = %.3f', pval1_2_window2, pval1_3_window2, pval1_4_window2, pval2_3_window2, pval2_4_window2, pval3_4_window2), 'FontSize', 12);
    
    % Display the actual values for verification
    disp(['Panel m - Normalized Mean Z-Score difference for -4 to 0 s: Fast = ', num2str(meanZScoresFastestWindow1_m), ', Slow = ', num2str(meanZScoresSlowestWindow1_m), ', Near Miss = ', num2str(meanZScoresFastestMissWindow1_m), ', False Start = ', num2str(meanZScoresSlowestMissWindow1_m)]);
    disp(['Panel m - Mean Z-Score difference for 0 to 2 s: Fast = ', num2str(meanZScoresFastestWindow2_m), ', Slow = ', num2str(meanZScoresSlowestWindow2_m), ', Near Miss = ', num2str(meanZScoresFastestMissWindow2_m), ', False Start = ', num2str(meanZScoresSlowestMissWindow2_m)]);


    %SUPPLEMENTAL FIGURE 3 BAR PLOTS
    % Define RT ranges for Fast and Slow Trials
    fastestStart = 0.02;
    fastestEnd = 0.03;
    slowestStart = median(sortedRT);
    slowestEnd = 0.5;

    % Identify the indices of trials within the defined RT ranges for fastest and slowest trials
    fastestIndices = (combinedRTToneSuccessful >= fastestStart) & (combinedRTToneSuccessful < fastestEnd);
    slowestIndices = (combinedRTToneSuccessful >= slowestStart) & (combinedRTToneSuccessful <= slowestEnd);

    % Define time windows for analysis
    timeWindowPreparatory = (psthTime >= -4) & (psthTime < 0);
    timeWindowAnticipatory = (psthTime >= 0) & (psthTime <= 2);
    
    % Initialize arrays for storing data
    diffs = zeros(1, 4);
    sems = zeros(1, 4);
    pvals = zeros(1, 4);
    fiberData = cell(1, 2);
    
    % Loop through both fibers
    for fiberIndex = 1:2
        % Extract z-score data for the fastest and slowest trials for current fiber
        zScoresFastest = combinedPsthSuccessMatrix(fastestIndices, :, fiberIndex);
        zScoresSlowest = combinedPsthSuccessMatrix(slowestIndices, :, fiberIndex);
    
        % Extract z-scores for the defined time windows
        zScoresFastestPrep = zScoresFastest(:, timeWindowPreparatory);
        zScoresSlowestPrep = zScoresSlowest(:, timeWindowPreparatory);
        
        zScoresFastestAntic = zScoresFastest(:, timeWindowAnticipatory);
        zScoresSlowestAntic = zScoresSlowest(:, timeWindowAnticipatory);
        
        % Compute the mean z-scores for the time windows
        meanZScoresFastestPrep = mean(zScoresFastestPrep(:));
        meanZScoresSlowestPrep = mean(zScoresSlowestPrep(:));
        
        meanZScoresFastestAntic = mean(zScoresFastestAntic(:));
        meanZScoresSlowestAntic = mean(zScoresSlowestAntic(:));
        
        % Calculate the differences in z-scores between fast and slow trials for each time window
        diffPrep = meanZScoresFastestPrep - meanZScoresSlowestPrep;
        diffAntic = meanZScoresFastestAntic - meanZScoresSlowestAntic;
        
        % Store data for later comparisons
        diffs(2 * (fiberIndex - 1) + 1) = diffPrep;
        diffs(2 * (fiberIndex - 1) + 2) = diffAntic;
        
        % Calculate SEMs
        semPrep = std(zScoresFastestPrep(:)) / sqrt(length(zScoresFastestPrep(:))) + ...
                  std(zScoresSlowestPrep(:)) / sqrt(length(zScoresSlowestPrep(:)));
        semAntic = std(zScoresFastestAntic(:)) / sqrt(length(zScoresFastestAntic(:))) + ...
                   std(zScoresSlowestAntic(:)) / sqrt(length(zScoresSlowestAntic(:)));
        
        sems(2 * (fiberIndex - 1) + 1) = semPrep;
        sems(2 * (fiberIndex - 1) + 2) = semAntic;
        
        % Perform two-sample t-tests
        [~, pvals(2 * (fiberIndex - 1) + 1)] = ttest2(zScoresFastestPrep(:), zScoresSlowestPrep(:));
        [~, pvals(2 * (fiberIndex - 1) + 2)] = ttest2(zScoresFastestAntic(:), zScoresSlowestAntic(:));
        
        % Store z-score data for fiber comparisons
        fiberData{fiberIndex} = {zScoresFastestPrep, zScoresSlowestPrep, zScoresFastestAntic, zScoresSlowestAntic};
    end
    
    % Extract data for fiber 1 and fiber 2
    fiber1_fast_preparatory = fiberData{1, 1}{1};  % 87x80 double
    fiber1_slow_preparatory = fiberData{1, 1}{2};  % 1199x80 double
    fiber1_fast_anticipatory = fiberData{1, 1}{3}; % 87x41 double
    fiber1_slow_anticipatory = fiberData{1, 1}{4}; % 1199x41 double
    
    fiber2_fast_preparatory = fiberData{1, 2}{1};  % 87x80 double
    fiber2_slow_preparatory = fiberData{1, 2}{2};  % 1199x80 double
    fiber2_fast_anticipatory = fiberData{1, 2}{3}; % 87x41 double
    fiber2_slow_anticipatory = fiberData{1, 2}{4}; % 1199x41 double
    
    % Compute the mean z-score of slow trials during preparatory period for each fiber
    mean_slow_preparatory_fiber1 = mean(fiber1_slow_preparatory, 1);
    mean_slow_preparatory_fiber2 = mean(fiber2_slow_preparatory, 1);
    
    % Subtract the mean of slow trials from each individual fast trial for preparatory period
    diff_preparatory_fiber1 = fiber1_fast_preparatory - mean_slow_preparatory_fiber1;
    diff_preparatory_fiber2 = fiber2_fast_preparatory - mean_slow_preparatory_fiber2;
    
    % Perform t-test for preparatory period
    [h_prep, p_prep] = ttest2(diff_preparatory_fiber1(:), diff_preparatory_fiber2(:));
    
    % Compute the mean z-score of slow trials during anticipatory period for each fiber
    mean_slow_anticipatory_fiber1 = mean(fiber1_slow_anticipatory, 1);
    mean_slow_anticipatory_fiber2 = mean(fiber2_slow_anticipatory, 1);
    
    % Subtract the mean of slow trials from each individual fast trial for anticipatory period
    diff_anticipatory_fiber1 = fiber1_fast_anticipatory - mean_slow_anticipatory_fiber1;
    diff_anticipatory_fiber2 = fiber2_fast_anticipatory - mean_slow_anticipatory_fiber2;
    
    % Perform t-test for anticipatory period
    [h_anticipatory, p_anticipatory] = ttest2(diff_anticipatory_fiber1(:), diff_anticipatory_fiber2(:));
    
    % Create a single plot with all four bars and set dimensions for a square plot
    figure('Position', [100, 100, 700, 700]);  % Define a square dimension 
    barData = [diffs(1), diffs(3), diffs(2), diffs(4)];
    bar(1:4, barData, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on;
    
    % Add error bars (SEM)
    errorbar(1:4, barData, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    % Annotate p-values for each bar comparison
    text(1, barData(1) + 0.02, sprintf('p=%.3f', pvals(1)), 'HorizontalAlignment', 'center');
    text(2, barData(2) + 0.02, sprintf('p=%.3f', pvals(3)), 'HorizontalAlignment', 'center');
    text(3, barData(3) + 0.02, sprintf('p=%.3f', pvals(2)), 'HorizontalAlignment', 'center');
    text(4, barData(4) + 0.02, sprintf('p=%.3f', pvals(4)), 'HorizontalAlignment', 'center');
    
    % Customize plot appearance
    set(gca, 'XTick', 1:4, 'XTickLabel', {'MD → dmPFC Preparatory Period (-4s to 0s)', 'MD Preparatory Period (-4s to 0s)', 'MD → dmPFC Anticipatory Period (0s to 2s)', 'MD Anticipatory Period (0s to 2s)'});
    ylabel('Difference in Mean Z-Score');
    title('Difference in Z-Scores Between Fast and Slow Trials for Each Fiber');
    ylim([-0.25 0.25]);
    box off;
    grid off;
    
    % Display p-values comparing fibers and periods
    annotation('textbox', [0.1, 0.9, 0.8, 0.1], 'String', ...
               sprintf('MD → dmPFC vs MD Preparatory p=%.3f | MD → dmPFC vs MD Anticipatory p=%.3f', ...
                       p_prep, p_anticipatory), ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end
