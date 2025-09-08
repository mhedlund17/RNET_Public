%% MFP Analysis + Manuscript Figures
% Date: 2025-09-05
% Requires: processDataForExperiment_zscored.m

%% -------------------- Parameters --------------------
% Clean the workspace and close all figures
clear;
close all;

% change basePath based on location of data
basePath = '/Users/sina/Documents/Stanford/Misc/D-Lab/mfp_analysis/data/';
protocols = {'test'};
mice = {'m1','m2', 'm3', 'm4', 'm5', 'm6'};
first_experiment_count = 0;
last_experiment_count  = 34;
numfibers = 2;          % num photometry inputs, keep at 2

%% -------------------- Load/Aggregate --------------------

% Initialize combined variables
combinedpsthTime_struct = struct();
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
combinedpsthTime = [];
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

        combinedpsthTime_struct.(mice{mouseIdx}) = [];
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
            combinedpsthTime_struct.(mice{mouseIdx}) = [combinedpsthTime_struct.(mice{mouseIdx}); psthTime];
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
        mouseData_combinedpsthTime = combinedpsthTime_struct.(mouseName);
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
        combinedpsthTime = cat(1, combinedpsthTime, mouseData_combinedpsthTime);
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

psthTime = combinedpsthTime(1,:);

%anticipatory and preparatory definition
timeWindow1 = (psthTime >= -4) & (psthTime < 0);
timeWindow2 = (psthTime >= 0)  & (psthTime <= 2);

%% =========================================================
%% ===================== FIGURE 4 ==========================
%% =========================================================

%% =================== Figure 4C ===================
figure('Name','Figure 4C — Distribution of Successful RTs');
histogram(combinedRTToneSuccessful, 'BinWidth', 0.01);
xlabel('Reaction Time (s)'); ylabel('Count'); title('Successful RTs'); box off;

%% =================== Figure 4D ===================
% 2nd Fastest vs 2nd Slowest (Fiber 1)
sortedRTs = sort(combinedRTToneSuccessful);
if numel(sortedRTs) >= 2
    rt2fast = sortedRTs(2);
    rt2slow = sortedRTs(end-1);
    idx2fast = find(combinedRTToneSuccessful == rt2fast, 1, 'first');
    idx2slow = find(combinedRTToneSuccessful == rt2slow, 1, 'first');

    yFast = squeeze(combinedPsthSuccessMatrix(idx2fast, :, 1));
    ySlow = squeeze(combinedPsthSuccessMatrix(idx2slow, :, 1));

    figure('Name','Figure 4D — 2nd Fastest (Blue) vs 2nd Slowest (Red)');
    plot(psthTime, yFast, 'b', 'LineWidth', 1.5); hold on;
    plot(psthTime, ySlow, 'r', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Z-Score'); title('Example Trials (Fiber 1)');
    xlim([-10, 10]); ylim([-2, 2.5]);
    plot([0 0], ylim, 'k-'); plot([2 2], ylim, 'k--'); plot([2.5 2.5], ylim, 'k--');
    legend({'2nd Fastest','2nd Slowest'}, 'Location','best'); box off; hold off;

    disp(['[Fig 4D] 2nd Fastest RT: ', num2str(rt2fast)]);
    disp(['[Fig 4D] 2nd Slowest RT: ', num2str(rt2slow)]);
end

%% -------------------- Precompute adaptive RT bin edges (for Fig 4E–F) --------------------
% Define minimum number of observations per bin
minObservations = 100;

% Initialize variables
binStart = 0;
binEnd = 0.02; % 20 ms bins
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

%% =================== Figure 4E ===================
% Adaptive RT-Binned Gradient PSTH (Successes)
figure('Name','Figure 4E — Adaptive RT-Binned Gradient PSTH (Successes)');
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

%% =================== Figure 4F ===================    
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
[yfit, ~] = polyval(p, RTBins, S);
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

% Calculate p value
SXX      = sum( (RTBins - mean(RTBins)).^2 );   % variance of predictor
df       = max(1, n-2);                         % degrees of freedom
s2       = SSresid / df;                        % residual variance
se_slope = sqrt( s2 / max(eps, SXX) );          % SE of slope
t_slope  = p(1) / se_slope;                     % t-stat for slope (p(1) is slope from polyfit)
pval     = 2 * tcdf(-abs(t_slope), df);         % two-sided p-value

if pval < 0.001
    pval_text = 'p < 0.001';
else
    pval_text = sprintf('p = %.3f', pval);
end

% Add R^2 and p-value to the plot
text(0.3, max(meanZScores) + 0.05, sprintf('R^2 = %.2f\n%s', rsq, pval_text), 'FontSize', 12, 'HorizontalAlignment', 'center');

xlabel('Mean RT (s)');
ylabel('Mean z-score from 1 to 2s');
title('Scatter plot of RT bins vs. mean z-score from 1 to 2s');

% Adjust the axis limits if necessary
xlim([0, max(RTBins) + 0.01]);
ylim([min(meanZScores) - 0.05, max(meanZScores) + 0.05]);

grid off;
hold off;
box off;

%% =================== Figure 4G ===================
% Gradient PSTH by Fastest Segments (1/10 → 1/2)
segments = [10 9 8 7 6 5 4 3 2];
sortedRT = sort(combinedRTToneSuccessful, 'ascend');
N = numel(sortedRT);
colorsG = jet(numel(segments));

figure('Name','Figure 4G — Gradient PSTH by Fastest Segments');
for p = 1:numfibers
    subplot(numfibers,1,p); hold on;
    for s = 1:numel(segments)
        idxEnd = max(1, floor(N/segments(s)));
        thr = sortedRT(idxEnd);
        segMask = combinedRTToneSuccessful <= thr;

        psthCurrentSegment = combinedPsthSuccessMatrix(segMask, :, p);
        if isempty(psthCurrentSegment), continue; end
        meanCurrentSegment = mean(psthCurrentSegment, 1);
        errCurrentSegment  = std(psthCurrentSegment, 0, 1) / sqrt(size(psthCurrentSegment, 1));

        plot(psthTime, meanCurrentSegment, 'color', colorsG(s,:), 'LineWidth', 1.5);
        fill([psthTime, fliplr(psthTime)], ...
             [meanCurrentSegment + errCurrentSegment, fliplr(meanCurrentSegment - errCurrentSegment)], ...
             colorsG(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    xlabel('Time (s)'); ylabel('Response (z)');
    title(['Fastest Segments (Fiber ', num2str(p), ')']);
    plot([0 0], [-1 2.5], 'k-'); plot([2 2], [-1 2.5], 'k--'); plot([2.5 2.5], [-1 2.5], 'k--');
    ylim([-0.75 1.75]); box off; hold off;
end

%% =================== Figure 4H ===================
% Calculate the indices for each segment
numTrials = length(sortedRT);
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

grid off; hold off; box off;

%% =================== Figure 4I ===================
% Define variables
fastestStart = 0;       % default "fast" lower bound
fastestEnd   = 0.02;    % default "fast" upper bound (0–20 ms)
slowestStart = fastestEnd;
slowestEnd   = 0.5;

% Overlay PSTH: Fast (<=20 ms) vs Slow (>20–500 ms), Successes
fastMask = (combinedRTToneSuccessful >= fastestStart) & (combinedRTToneSuccessful < fastestEnd);
slowMask = (combinedRTToneSuccessful >= slowestStart) & (combinedRTToneSuccessful <= slowestEnd);

figure('Name','Figure 4I — Overlay PSTH: Fast vs Slow (Successes)');
for p = 1:numfibers
    subplot(numfibers,1,p);
    psthFast = combinedPsthSuccessMatrix(fastMask, :, p);
    psthSlow = combinedPsthSuccessMatrix(slowMask, :, p);

    mFast = mean(psthFast, 1);
    eFast = std(psthFast, 0, 1) / max(1, sqrt(size(psthFast,1)));
    mSlow = mean(psthSlow, 1);
    eSlow = std(psthSlow, 0, 1) / max(1, sqrt(size(psthSlow,1)));

    plot(psthTime, mFast, 'b'); hold on;
    if size(psthFast,1) > 1
        fill([psthTime, fliplr(psthTime)], [mFast+eFast, fliplr(mFast-eFast)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    plot(psthTime, mSlow, 'r');
    if size(psthSlow,1) > 1
        fill([psthTime, fliplr(psthTime)], [mSlow+eSlow, fliplr(mSlow-eSlow)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    xlabel('Time (s)'); ylabel('Response (z)');
    title(['Fast (N=', num2str(size(psthFast,1)), ', meanRT=', num2str(mean(combinedRTToneSuccessful(fastMask))), ...
           ' s) vs Slow (N=', num2str(size(psthSlow,1)), ', meanRT=', num2str(mean(combinedRTToneSuccessful(slowMask))), ...
           ' s) — Fiber ', num2str(p)]);
    plot([0 0], [-1 2.5], 'k-'); plot([2 2], [-1 2.5], 'k--'); plot([2.5 2.5], [-1 2.5], 'k--');
    ylim([-0.75 1.75]); box off; hold off;
end

%% =================== Figure 4J ===================
% Δ Mean z (Fast − Slow), Two Windows (Successes, Fiber 1)
zFast_w1 = combinedPsthSuccessMatrix(fastMask, timeWindow1, 1);
zSlow_w1 = combinedPsthSuccessMatrix(slowMask, timeWindow1, 1);
zFast_w2 = combinedPsthSuccessMatrix(fastMask, timeWindow2, 1);
zSlow_w2 = combinedPsthSuccessMatrix(slowMask, timeWindow2, 1);

diff1 = mean(zFast_w1(:)) - mean(zSlow_w1(:));
diff2 = mean(zFast_w2(:)) - mean(zSlow_w2(:));
[~, pJ1] = ttest2(zFast_w1(:), zSlow_w1(:));
[~, pJ2] = ttest2(zFast_w2(:), zSlow_w2(:));

sem1 = hypot(std(zFast_w1(:))/sqrt(max(1,numel(zFast_w1))), std(zSlow_w1(:))/sqrt(max(1,numel(zSlow_w1))));
sem2 = hypot(std(zFast_w2(:))/sqrt(max(1,numel(zFast_w2))), std(zSlow_w2(:))/sqrt(max(1,numel(zSlow_w2))));

figure('Name','Figure 4J — Δ Mean z (Fast − Slow)');
bar([1 2], [diff1 diff2], 'FaceColor','none', 'EdgeColor','k', 'LineWidth',1.5); hold on;
errorbar([1 2], [diff1 diff2], [sem1 sem2], 'k', 'LineStyle','none', 'LineWidth',1.5);
set(gca,'XTick',[1 2],'XTickLabel',{'−4 to 0 s','0 to 2 s'});
ylabel('Δ Mean z (Fast − Slow)'); title('Fast vs Slow (Successes)');
text(1, diff1 + sem1 + 0.02, sprintf('p = %.3g', pJ1), 'HorizontalAlignment','center');
text(2, diff2 + sem2 + 0.02, sprintf('p = %.3g', pJ2), 'HorizontalAlignment','center');
box off; hold off;

%% =================== Figure 4K ===================
% Need to load 100ms Data for this panel
% Define variables
RTdivision   = 3; % divide the data into thirds
idxFastest   = floor(length(combinedRTToneSuccessful) / RTdivision);
idxSlowest   = floor(length(combinedRTToneSuccessful) - idxFastest);
fastestEnd   = sortedRT(idxFastest);
slowestStart = sortedRT(idxSlowest);
fastestStart = 0;
slowestEnd   = 0.1; 

% Overlay PSTH: Fast (<=20 ms) vs Slow (>20–500 ms), Successes
fastMask = (combinedRTToneSuccessful >= fastestStart) & (combinedRTToneSuccessful < fastestEnd);
slowMask = (combinedRTToneSuccessful >= slowestStart) & (combinedRTToneSuccessful <= slowestEnd);

figure('Name','Figure 4I — Overlay PSTH: Fast vs Slow (Successes)');
for p = 1:numfibers
    subplot(numfibers,1,p);
    psthFast = combinedPsthSuccessMatrix(fastMask, :, p);
    psthSlow = combinedPsthSuccessMatrix(slowMask, :, p);

    mFast = mean(psthFast, 1);
    eFast = std(psthFast, 0, 1) / max(1, sqrt(size(psthFast,1)));
    mSlow = mean(psthSlow, 1);
    eSlow = std(psthSlow, 0, 1) / max(1, sqrt(size(psthSlow,1)));

    plot(psthTime, mFast, 'b'); hold on;
    if size(psthFast,1) > 1
        fill([psthTime, fliplr(psthTime)], [mFast+eFast, fliplr(mFast-eFast)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    plot(psthTime, mSlow, 'c');
    if size(psthSlow,1) > 1
        fill([psthTime, fliplr(psthTime)], [mSlow+eSlow, fliplr(mSlow-eSlow)], 'c', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    xlabel('Time (s)'); ylabel('Response (z)');
    title(['Fast (N=', num2str(size(psthFast,1)), ', meanRT=', num2str(mean(combinedRTToneSuccessful(fastMask))), ...
           ' s) vs Slow (N=', num2str(size(psthSlow,1)), ', meanRT=', num2str(mean(combinedRTToneSuccessful(slowMask))), ...
           ' s) — Fiber ', num2str(p)]);
    plot([0 0], [-1 2.5], 'k-'); plot([2 2], [-1 2.5], 'k--'); plot([2.5 2.5], [-1 2.5], 'k--');
    ylim([-0.75 1.75]); box off; hold off;
end

%% =================== Figure 4L ===================
% Need to load 100ms Data for this panel
% Within 0–100 ms: Δ Mean z (Fastest Third vs Slowest Third)
sortedRT_f = sort(combinedRTToneSuccessful);
idxFastest_f = max(1, floor(numel(sortedRT_f)/3));
idxSlowest_f = max(1, floor(numel(sortedRT_f) - idxFastest_f));
fastestEnd_f = sortedRT_f(idxFastest_f);
slowestStart_f = sortedRT_f(idxSlowest_f);
slowestEnd_f = 0.1;

fastestIndices_f = (combinedRTToneSuccessful >= 0) & (combinedRTToneSuccessful <  fastestEnd_f);
slowestIndices_f = (combinedRTToneSuccessful >= slowestStart_f) & (combinedRTToneSuccessful <= slowestEnd_f);

zF1 = combinedPsthSuccessMatrix(fastestIndices_f, timeWindow1, 1);
zS1 = combinedPsthSuccessMatrix(slowestIndices_f, timeWindow1, 1);
zF2 = combinedPsthSuccessMatrix(fastestIndices_f, timeWindow2, 1);
zS2 = combinedPsthSuccessMatrix(slowestIndices_f, timeWindow2, 1);

d1 = mean(zF1(:)) - mean(zS1(:));
d2 = mean(zF2(:)) - mean(zS2(:));
[~, pL1] = ttest2(zF1(:), zS1(:));
[~, pL2] = ttest2(zF2(:), zS2(:));

semL1 = hypot(std(zF1(:))/sqrt(max(1,numel(zF1))), std(zS1(:))/sqrt(max(1,numel(zS1))));
semL2 = hypot(std(zF2(:))/sqrt(max(1,numel(zF2))), std(zS2(:))/sqrt(max(1,numel(zS2))));

figure('Name','Figure 4L — Δ Mean z (0–100 ms subset)');
bar([1 2], [d1 d2], 'FaceColor','none','EdgeColor','k','LineWidth',1.5); hold on;
errorbar([1 2], [d1 d2], [semL1 semL2], 'k','LineStyle','none','LineWidth',1.5);
set(gca,'XTick',[1 2],'XTickLabel',{'−4 to 0 s','0 to 2 s'});
ylabel('Δ Mean z (Fastest Third − Slowest Third)'); title('0–100 ms Subset');
text(1, d1 + semL1 + 0.02, sprintf('p = %.3g', pL1), 'HorizontalAlignment','center');
text(2, d2 + semL2 + 0.02, sprintf('p = %.3g', pL2), 'HorizontalAlignment','center');
box off; hold off;

%% =========================================================
%% ================Supplemental Figure S7===================
%% =========================================================
% Adjust fastestEnd based on desired time window
% Define variables
fastestStart = 0;      
fastestEnd   = 0.1;    
slowestStart = fastestEnd;
slowestEnd   = 0.5;

% Overlay PSTH: Fast (<=20 ms) vs Slow (>20–500 ms), Successes
fastMask = (combinedRTToneSuccessful >= fastestStart) & (combinedRTToneSuccessful < fastestEnd);
slowMask = (combinedRTToneSuccessful >= slowestStart) & (combinedRTToneSuccessful <= slowestEnd);

figure('Name','Figure S7 — Overlay PSTH: Fast vs Slow (Successes)');
for p = 1:numfibers
    subplot(numfibers,1,p);
    psthFast = combinedPsthSuccessMatrix(fastMask, :, p);
    psthSlow = combinedPsthSuccessMatrix(slowMask, :, p);

    mFast = mean(psthFast, 1);
    eFast = std(psthFast, 0, 1) / max(1, sqrt(size(psthFast,1)));
    mSlow = mean(psthSlow, 1);
    eSlow = std(psthSlow, 0, 1) / max(1, sqrt(size(psthSlow,1)));

    plot(psthTime, mFast, 'b'); hold on;
    if size(psthFast,1) > 1
        fill([psthTime, fliplr(psthTime)], [mFast+eFast, fliplr(mFast-eFast)], 'b', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    plot(psthTime, mSlow, 'r');
    if size(psthSlow,1) > 1
        fill([psthTime, fliplr(psthTime)], [mSlow+eSlow, fliplr(mSlow-eSlow)], 'r', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    xlabel('Time (s)'); ylabel('Response (z)');
    title(['Fast (N=', num2str(size(psthFast,1)), ', meanRT=', num2str(mean(combinedRTToneSuccessful(fastMask))), ...
         ' s) vs Slow (N=', num2str(size(psthSlow,1)), ', meanRT=', num2str(mean(combinedRTToneSuccessful(slowMask))), ...
          ' s) — Fiber ', num2str(p)]);
    plot([0 0], [-1 2.5], 'k-'); plot([2 2], [-1 2.5], 'k--'); plot([2.5 2.5], [-1 2.5], 'k--');
    ylim([-0.75 1.75]); box off; hold off;
end

% Barplots
% Δ Mean z (Fast − Slow), Two Windows (Successes, Fiber 1)
zFast_w1 = combinedPsthSuccessMatrix(fastMask, timeWindow1, 1);
zSlow_w1 = combinedPsthSuccessMatrix(slowMask, timeWindow1, 1);
zFast_w2 = combinedPsthSuccessMatrix(fastMask, timeWindow2, 1);
zSlow_w2 = combinedPsthSuccessMatrix(slowMask, timeWindow2, 1);

diff1 = mean(zFast_w1(:)) - mean(zSlow_w1(:));
diff2 = mean(zFast_w2(:)) - mean(zSlow_w2(:));
[~, pJ1] = ttest2(zFast_w1(:), zSlow_w1(:));
[~, pJ2] = ttest2(zFast_w2(:), zSlow_w2(:));

sem1 = hypot(std(zFast_w1(:))/sqrt(max(1,numel(zFast_w1))), std(zSlow_w1(:))/sqrt(max(1,numel(zSlow_w1))));
sem2 = hypot(std(zFast_w2(:))/sqrt(max(1,numel(zFast_w2))), std(zSlow_w2(:))/sqrt(max(1,numel(zSlow_w2))));

figure('Name','Figure S7 — Δ Mean z (Fast − Slow)');
bar([1 2], [diff1 diff2], 'FaceColor','none', 'EdgeColor','k', 'LineWidth',1.5); hold on;
errorbar([1 2], [diff1 diff2], [sem1 sem2], 'k', 'LineStyle','none', 'LineWidth',1.5);
set(gca,'XTick',[1 2],'XTickLabel',{'−4 to 0 s','0 to 2 s'});
ylabel('Δ Mean z (Fast − Slow)'); title('Fast vs Slow (Successes)');
text(1, diff1 + sem1 + 0.02, sprintf('p = %.3g', pJ1), 'HorizontalAlignment','center');
text(2, diff2 + sem2 + 0.02, sprintf('p = %.3g', pJ2), 'HorizontalAlignment','center');
ylim([-0.1 0.2]); box off; hold off;

%% =========================================================
%% ================Supplemental Figure S8===================
%% =========================================================
% Fast vs Slow ≥ median, by Fiber & Window
% Adjust fastestStart_S8 and slowestStart_S8 to get desired time window
fastestStart_S8 = 0.0; fastestEnd_S8 = 0.01;
slowestStart_S8 = median(combinedRTToneSuccessful); slowestEnd_S8 = slowestEnd;

idxFast_S8 = (combinedRTToneSuccessful >= fastestStart_S8) & (combinedRTToneSuccessful < fastestEnd_S8);
idxSlow_S8 = (combinedRTToneSuccessful >= slowestStart_S8) & (combinedRTToneSuccessful <= slowestEnd_S8);

tw1 = timeWindow1; tw2 = timeWindow2; 

diffs = zeros(1,4); sems = zeros(1,4); pvals = zeros(1,4);
for p = 1:min(2, numfibers)
    zF_prep = combinedPsthSuccessMatrix(idxFast_S8, tw1, p);
    zS_prep = combinedPsthSuccessMatrix(idxSlow_S8, tw1, p);
    zF_anti = combinedPsthSuccessMatrix(idxFast_S8, tw2, p);
    zS_anti = combinedPsthSuccessMatrix(idxSlow_S8, tw2, p);

    diffs(2*(p-1)+1) = mean(zF_prep(:)) - mean(zS_prep(:));
    diffs(2*(p-1)+2) = mean(zF_anti(:)) - mean(zS_anti(:));

    sems(2*(p-1)+1) = hypot(std(zF_prep(:))/sqrt(max(1,numel(zF_prep))), std(zS_prep(:))/sqrt(max(1,numel(zS_prep))));
    sems(2*(p-1)+2) = hypot(std(zF_anti(:))/sqrt(max(1,numel(zF_anti))), std(zS_anti(:))/sqrt(max(1,numel(zS_anti))));

    [~, pvals(2*(p-1)+1)] = ttest2(zF_prep(:), zS_prep(:));
    [~, pvals(2*(p-1)+2)] = ttest2(zF_anti(:), zS_anti(:));
end

if numfibers >= 2
    meanSlowPrep1 = mean(combinedPsthSuccessMatrix(idxSlow_S8, tw1, 1), 1);
    meanSlowPrep2 = mean(combinedPsthSuccessMatrix(idxSlow_S8, tw1, 2), 1);
    meanSlowAnti1 = mean(combinedPsthSuccessMatrix(idxSlow_S8, tw2, 1), 1);
    meanSlowAnti2 = mean(combinedPsthSuccessMatrix(idxSlow_S8, tw2, 2), 1);

    diffPrep1 = combinedPsthSuccessMatrix(idxFast_S8, tw1, 1) - meanSlowPrep1;
    diffPrep2 = combinedPsthSuccessMatrix(idxFast_S8, tw1, 2) - meanSlowPrep2;
    diffAnti1 = combinedPsthSuccessMatrix(idxFast_S8, tw2, 1) - meanSlowAnti1;
    diffAnti2 = combinedPsthSuccessMatrix(idxFast_S8, tw2, 2) - meanSlowAnti2;

    [~, p_prepFib] = ttest2(diffPrep1(:), diffPrep2(:));
    [~, p_antiFib] = ttest2(diffAnti1(:), diffAnti2(:));
else
    p_prepFib = NaN; p_antiFib = NaN;
end

figure('Name','Supplemental Figure S8 — Δ Mean z (Fast 20–30 ms − Slow ≥ median)');
bar(1:4, [diffs(1), diffs(3), diffs(2), diffs(4)], 'FaceColor','none','EdgeColor','k','LineWidth',1.5); hold on;
errorbar(1:4, [diffs(1), diffs(3), diffs(2), diffs(4)], [sems(1), sems(3), sems(2), sems(4)], 'k', 'LineStyle','none', 'LineWidth',1.5);
set(gca, 'XTick', 1:4, 'XTickLabel', {'MD→dmPFC (−4–0)', 'MD (−4–0)', 'MD→dmPFC (0–2)', 'MD (0–2)'});
ylabel('Δ Mean z (Fast − Slow)'); title('Fast 20–30 ms vs Slow ≥ median');
text(1, diffs(1)+sems(1)+0.02, sprintf('p=%.3g', pvals(1)), 'HorizontalAlignment','center');
text(2, diffs(3)+sems(3)+0.02, sprintf('p=%.3g', pvals(3)), 'HorizontalAlignment','center');
text(3, diffs(2)+sems(2)+0.02, sprintf('p=%.3g', pvals(2)), 'HorizontalAlignment','center');
text(4, diffs(4)+sems(4)+0.02, sprintf('p=%.3g', pvals(4)), 'HorizontalAlignment','center');
ylim([-0.25 0.25]); box off; hold off;

annotation('textbox',[0.1 0.9 0.8 0.1],'String', ...
    sprintf('Cross-fiber: prep p=%.3g | anti p=%.3g', p_prepFib, p_antiFib), ...
    'EdgeColor','none','HorizontalAlignment','center');

%% =========================================================
%% ================Supplemental Figure S9===================
%% =========================================================
%% ===============Panel S9A===================
% Define Variables
fastestStart = 0;       % default "fast" lower bound
fastestEnd   = 0.02;    % default "fast" upper bound (0–20 ms)
slowestStart = fastestEnd;
slowestEnd   = 0.5;
NearMissThreshold = 1.98; % near misses (green), seconds
FalseStartThreshold = 0.1; % false starts (magenta), seconds

% Indices for the fastest, slowest, near miss, and false start
fastestIndices = combinedRTToneSuccessful >= fastestStart & combinedRTToneSuccessful < fastestEnd;
slowestIndices = combinedRTToneSuccessful > slowestStart & combinedRTToneSuccessful <= slowestEnd;
nearMissIndices = combinedRTToneUnsuccessful > NearMissThreshold;
falseStartIndices = combinedRTToneUnsuccessful < FalseStartThreshold;

% PSTH for the fastest, slowest, near miss, and false start
psthFastest = combinedPsthSuccessMatrix(fastestIndices, :, :);
psthSlowest = combinedPsthSuccessMatrix(slowestIndices, :, :);
psthNearMiss = combinedPsthWithoutSecondMatrix(nearMissIndices, :, :);
psthFalseStarts = combinedPsthWithoutSecondMatrix(falseStartIndices, :, :);
            
figure; 
for p = 1:numfibers
    subplot(numfibers, 1, p);
    
    % Compute the mean and error for fastest data
    meanFastest = mean(psthFastest(:,:,p), 1);
    errFastest = std(psthFastest(:,:,p), 0, 1) / sqrt(size(psthFastest, 1));
    
    % Compute the mean and error for slowest data
    meanSlowest = mean(psthSlowest(:,:,p),1);
    errSlowest = std(psthSlowest(:,:,p),0,1) / sqrt(size(psthSlowest, 1));

    % Compute the mean and error for near miss trials
    meanNearMiss = mean(psthNearMiss(:,:,p),1);
    errNearMiss = std(psthNearMiss(:,:,p),0,1) / sqrt(size(psthNearMiss, 1));
    
    % Compute the mean and error for false start trials
    meanFalseStart = mean(psthFalseStarts(:,:,p),1);
    errFalseStart = std(psthFalseStarts(:,:,p),0,1) / sqrt(size(psthFalseStarts, 1));

    % Plotting fastest trials
    plot(psthTime, meanFastest, 'b');
    hold on;
        fill([psthTime, fliplr(psthTime)], [meanFastest + errFastest, fliplr(meanFastest - errFastest)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    
    % Plotting slowest trials
    plot(psthTime, meanSlowest, 'r');
    if size(psthSlowest, 1) > 1 % Only plot error bands if there's more than one trial
        fill([psthTime, fliplr(psthTime)], [meanSlowest + errSlowest, fliplr(meanSlowest - errSlowest)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    
    % Plotting fastest starts 
    plot(psthTime, meanNearMiss, 'g');
    hold on;
    if size(psthNearMiss, 1) > 1 
        fill([psthTime, fliplr(psthTime)], [meanNearMiss + errNearMiss, fliplr(meanNearMiss - errNearMiss)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % Plotting slowest missed data
    plot(psthTime, meanFalseStart, 'm');
    if size(psthFalseStarts, 1) > 1 
        fill([psthTime, fliplr(psthTime)], [meanFalseStart + errFalseStart, fliplr(meanFalseStart - errFalseStart)], 'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    
    xlabel('Time (s)');
    ylabel('Response');
    title(['Overlay PSTH - Fastest RT Thirds (Blue), Slowest RT Thirds (Red), Near Misses (Green), False Starts (Magenta) - Fiber ', num2str(p)]);
    plot([0 0], [-1 2.5], 'k-', 'LineWidth', 2);
    plot([2 2], [-1 2.5], 'k--', 'LineWidth', 2);
    plot([2.5 2.5], [-1 2.5], 'k--', 'LineWidth', 2);
    box off;
    ylim([-0.75 1.75]);  % Set y-axis limits here
    hold off;
end

%% ===============Panel S9B===================
% Extract z-scores for the defined time windows
% Normalize the preparatory window values by adding 1
zScoresFastestWindow1_m = psthFastest(:, timeWindow1, 1) + 1;
zScoresSlowestWindow1_m = psthSlowest(:, timeWindow1, 1) + 1;
zScoresNearMissWindow1_m = psthNearMiss(:, timeWindow1, 1) + 1;
zScoresFalseStartWindow1_m = psthFalseStarts(:, timeWindow1, 1) + 1;

zScoresFastestWindow2_m = psthFastest(:, timeWindow2, 1);
zScoresSlowestWindow2_m = psthSlowest(:, timeWindow2, 1);
zScoresNearMissWindow2_m = psthNearMiss(:, timeWindow2, 1);
zScoresFalseStartWindow2_m = psthFalseStarts(:, timeWindow2, 1);

% Compute the mean z-scores for the time windows
meanZScoresFastestWindow1_m = mean(zScoresFastestWindow1_m(:));
meanZScoresSlowestWindow1_m = mean(zScoresSlowestWindow1_m(:));
meanZScoresNearMissWindow1_m = mean(zScoresNearMissWindow1_m(:));
meanZScoresFalseStartWindow1_m = mean(zScoresFalseStartWindow1_m(:));

meanZScoresFastestWindow2_m = mean(zScoresFastestWindow2_m(:));
meanZScoresSlowestWindow2_m = mean(zScoresSlowestWindow2_m(:));
meanZScoresNearMissWindow2_m = mean(zScoresNearMissWindow2_m(:));
meanZScoresFalseStartWindow2_m = mean(zScoresFalseStartWindow2_m(:));

% Calculate the SEM for the time windows
semFastestWindow1_m = std(zScoresFastestWindow1_m(:)) / sqrt(length(zScoresFastestWindow1_m(:)));
semSlowestWindow1_m = std(zScoresSlowestWindow1_m(:)) / sqrt(length(zScoresSlowestWindow1_m(:)));
semNearMissWindow1_m = std(zScoresNearMissWindow1_m(:)) / sqrt(length(zScoresNearMissWindow1_m(:)));
semFalseStartWindow1_m = std(zScoresFalseStartWindow1_m(:)) / sqrt(length(zScoresFalseStartWindow1_m(:)));

semFastestWindow2_m = std(zScoresFastestWindow2_m(:)) / sqrt(length(zScoresFastestWindow2_m(:)));
semSlowestWindow2_m = std(zScoresSlowestWindow2_m(:)) / sqrt(length(zScoresSlowestWindow2_m(:)));
semNearMissWindow2_m = std(zScoresFastestWindow2_m(:)) / sqrt(length(zScoresFastestWindow2_m(:)));
semFalseStartWindow2_m = std(zScoresSlowestWindow2_m(:)) / sqrt(length(zScoresSlowestWindow2_m(:)));

% Perform pairwise t-tests for time window (-4 to 0)
[~, pval1_2_window1] = ttest2(zScoresFastestWindow1_m(:), zScoresSlowestWindow1_m(:));
[~, pval1_3_window1] = ttest2(zScoresFastestWindow1_m(:), zScoresNearMissWindow1_m(:));
[~, pval1_4_window1] = ttest2(zScoresFastestWindow1_m(:), zScoresFalseStartWindow1_m(:));
[~, pval2_3_window1] = ttest2(zScoresSlowestWindow1_m(:), zScoresNearMissWindow1_m(:));
[~, pval2_4_window1] = ttest2(zScoresSlowestWindow1_m(:), zScoresFalseStartWindow1_m(:));
[~, pval3_4_window1] = ttest2(zScoresNearMissWindow1_m(:), zScoresFalseStartWindow1_m(:));

% Perform pairwise t-tests for time window (0 to 2)
[~, pval1_2_window2] = ttest2(zScoresFastestWindow2_m(:), zScoresSlowestWindow2_m(:));
[~, pval1_3_window2] = ttest2(zScoresFastestWindow2_m(:), zScoresNearMissWindow2_m(:));
[~, pval1_4_window2] = ttest2(zScoresFastestWindow2_m(:), zScoresFalseStartWindow2_m(:));
[~, pval2_3_window2] = ttest2(zScoresSlowestWindow2_m(:), zScoresNearMissWindow2_m(:));
[~, pval2_4_window2] = ttest2(zScoresSlowestWindow2_m(:), zScoresFalseStartWindow2_m(:));
[~, pval3_4_window2] = ttest2(zScoresNearMissWindow2_m(:), zScoresFalseStartWindow2_m(:));

% Create the bar plot for time window (-4 to 0)
figure;
subplot(1, 2, 1);
bar(1:4, [meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresNearMissWindow1_m, meanZScoresFalseStartWindow1_m], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on;
errorbar(1:4, [meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresNearMissWindow1_m, meanZScoresFalseStartWindow1_m], ...
         [semFastestWindow1_m, semSlowestWindow1_m, semNearMissWindow1_m, semFalseStartWindow1_m], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize bar colors
b = bar(1:4, [meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresNearMissWindow1_m, meanZScoresFalseStartWindow1_m], 'FaceColor', 'flat');
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
text(1.5, median([meanZScoresFastestWindow1_m, meanZScoresSlowestWindow1_m, meanZScoresNearMissWindow1_m, meanZScoresFalseStartWindow1_m]) + 0.1, sprintf('p_{1,2} = %.3f\np_{1,3} = %.3f\np_{1,4} = %.3f\np_{2,3} = %.3f\np_{2,4} = %.3f\np_{3,4} = %.3f', pval1_2_window1, pval1_3_window1, pval1_4_window1, pval2_3_window1, pval2_4_window1, pval3_4_window1), 'FontSize', 12);

% Create the bar plot for time window (0 to 2)
subplot(1, 2, 2);
bar(1:4, [meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresNearMissWindow2_m, meanZScoresFalseStartWindow2_m], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on;
errorbar(1:4, [meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresNearMissWindow2_m, meanZScoresFalseStartWindow2_m], ...
         [semFastestWindow2_m, semSlowestWindow2_m, semNearMissWindow2_m, semFalseStartWindow2_m], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Customize bar colors
b = bar(1:4, [meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresNearMissWindow2_m, meanZScoresFalseStartWindow2_m], 'FaceColor', 'flat');
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
text(1.5, median([meanZScoresFastestWindow2_m, meanZScoresSlowestWindow2_m, meanZScoresNearMissWindow2_m, meanZScoresFalseStartWindow2_m]) + 0.1, sprintf('p_{1,2} = %.3f\np_{1,3} = %.3f\np_{1,4} = %.3f\np_{2,3} = %.3f\np_{2,4} = %.3f\np_{3,4} = %.3f', pval1_2_window2, pval1_3_window2, pval1_4_window2, pval2_3_window2, pval2_4_window2, pval3_4_window2), 'FontSize', 12);

% Display the actual values
disp(['Panel m - Normalized Mean Z-Score difference for -4 to 0 s: Fast = ', num2str(meanZScoresFastestWindow1_m), ', Slow = ', num2str(meanZScoresSlowestWindow1_m), ', Near Miss = ', num2str(meanZScoresNearMissWindow1_m), ', False Start = ', num2str(meanZScoresFalseStartWindow1_m)]);
disp(['Panel m - Mean Z-Score difference for 0 to 2 s: Fast = ', num2str(meanZScoresFastestWindow2_m), ', Slow = ', num2str(meanZScoresSlowestWindow2_m), ', Near Miss = ', num2str(meanZScoresNearMissWindow2_m), ', False Start = ', num2str(meanZScoresFalseStartWindow2_m)]);

% % Save figures to the output directory
% outputDir = '/Users/sina/Desktop/untitled folder';
% figHandles = findobj('Type', 'figure');
% for i = 1:length(figHandles)
%     figHandle = figHandles(i);
%     originalPosition = get(figHandle, 'Position');
%     newPosition = originalPosition;
%     newPosition(3) = originalPosition(3) / 2.31; % Width
%     newPosition(4) = originalPosition(4) / 2.375; % Height
%     set(figHandle, 'Position', newPosition);
%     figName = sprintf('Figure_%d.svg', figHandle.Number);
%     saveas(figHandle, fullfile(outputDir, figName));
%     set(figHandle, 'Position', originalPosition);
% end