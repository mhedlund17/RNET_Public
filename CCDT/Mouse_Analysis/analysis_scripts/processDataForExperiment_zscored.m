function [psthMatrix, ITISuccessful, ITIUnsuccessful, RTToneSuccessful, RTToneUnsuccessful, RTSuccessful, RTUnsuccessful, psthSuccessMatrix, psthWithoutSecondMatrix, psthTime] = processDataForExperiment_zscored(filename)

% Check if file exists
if ~isfile(strcat(filename, '_bpod.mat'))
    disp(['File ', strcat(filename, '_bpod.mat'), ' does not exist. Skipping.']);
    psthMatrix = [];
    ITISuccessful = [];
    ITIUnsuccessful = [];
    RTToneSuccessful = [];
    RTToneUnsuccessful = [];
    RTSuccessful = [];
    RTUnsuccessful = [];
    psthSuccessMatrix = [];
    psthWithoutSecondMatrix = [];
    psthTime = [];
    return;
end

% Load data
ttl_input = 7;
numfibers = 2;
RTThresholdmax = 0.5; %Max threshold for successful lick, in seconds - keep at 0.5
RTThresholdmin = 0; %Min threshold for successful lick, in seconds - keep at 0
load(filename);
M = csvread(strcat(filename, '_logAI.csv'));  % Load the AIlog with TTLs and time marks
AItime = M(:, 1);  % Time tracking column
start_trial = M(:, ttl_input);  % TTL column

% Identify behavioral experiment start time from the TTLs
Time0 = start_trial > 1.5;  
listTime0 = find(Time0);
TrialStartTime = AItime(listTime0);

% Load bpod data and adjust times
bpod = load(strcat(filename, '_bpod.mat'));
bpodTrialStartTimes = bpod.SessionData.TrialStartTimestamp;
TrialStartAdj =  bpodTrialStartTimes + TrialStartTime(1) - bpodTrialStartTimes(1); % Adjust trial start times

% Pre-process Photometry Signals
sig = sig';
ref = ref'; 
sig(:, end-5:end, :) = [];
ref(:, end-5:end, :) = [];
fs = 1/20;  % Photometry sampling rate per channel
time = [0:fs:(length(sig)-1)*fs];
tPrev = 200; %200
tPost = 200; %200
psthMatrix = NaN(numel(TrialStartAdj), (tPrev+tPost+1), numfibers);
psthTime = [-tPrev*fs : fs : tPost*fs];

% Fiber Data Processing
for a = 1:numfibers
    p = polyfit(ref(a,:), sig(a,:), 1);
    ref_scaled = ref(a,:)*p(1) + p(2);
    
    % Plot raw and scaled signals
    figure(a);
    subplot(2, 1, 1);
    plot(time, sig(a,:));
    hold on;
    plot(time, ref_scaled, 'r');
    legend('Signal', 'Scaled control');
    title(['Fiber ', num2str(a), ': Signal and scaled control']);
    
    % Normalize signal and plot
    sig_norm(a, :) = sig(a, :) - ref_scaled + sig(a, 1);
    subplot(2, 1, 2);
    plot(time, sig_norm(a,:));
    for marker = 1:numel(TrialStartAdj)
        hold on;
        v = TrialStartAdj(marker);
        plot([v v], ylim);
    end
    title(['Fiber ', num2str(a), ': Normalized signal']);
end

% Post-processing: Z-score Calculation
baselineAvg=mean(sig_norm,2);
baselineStd=std(sig_norm,0,2);
for a = 1:numfibers
    sig_normZ(a,:)=(sig_norm(a,:)-baselineAvg(a))/baselineStd(a);
    psthArray = NaN(numel(TrialStartAdj), (tPrev+tPost+1)); 
    for x = 1:numel(TrialStartAdj)
        [~, thisTime] = min(abs(time-TrialStartAdj(x)));
        if (thisTime-tPrev) > 0 && (thisTime+tPost) <= length(sig_normZ)
            psthArray(x, :) = sig_normZ(a, thisTime-tPrev:thisTime+tPost);
        end
    end
    psthMatrix(:, :, a) = psthArray;
end

% Extract RT, ITI, and First Tone
validTrials = ~isnan(TrialStartAdj);
numTrials = length(validTrials);
RTSuccessful = NaN(1, numTrials);
RTUnsuccessful = NaN(1, numTrials);
ITI = NaN(1, numTrials);
FirstTone = NaN(1, numel(numTrials));

for i = 2:numTrials %Start from 2 to remove the first trial
    if isfield(bpod.SessionData.RawEvents.Trial{1, i}.Events, 'GlobalTimer2_Start')
        firstToneTime = bpod.SessionData.RawEvents.Trial{1,i}.Events.GlobalTimer2_Start;
        ITI(i) = bpod.SessionData.RawEvents.Trial{1, i}.States.ITI(:, 2);
        FirstTone(i) = bpod.SessionData.RawEvents.Trial{1,i}.Events.GlobalTimer2_Start; %FOR TRAIN/TEST PROTOCOL
        noITIlicks = true;
        for j = 1:length(bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In) %This will ignore licks before 10s of the first tone. Otherwise, successful trials will be defined as no licks until the second tone (with no licks during the ITI or first tone)
            if bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In(j) < firstToneTime && bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In(j) >= (firstToneTime - 10)
                noITIlicks = false;
                firstLickTime = bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In(j);
                break; 
            elseif bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In(j) >= firstToneTime
                firstLickTime = bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In(j);
                break; 
            end
        end
        %firstLickTime = bpod.SessionData.RawEvents.Trial{1, i}.Events.Port1In(1); % comment if uncommenting the above loop
        if isfield(bpod.SessionData.RawEvents.Trial{1, i}.Events, 'GlobalTimer3_Start')
            secondToneTime = bpod.SessionData.RawEvents.Trial{1, i}.Events.GlobalTimer3_Start;
            if firstLickTime > secondToneTime && noITIlicks
                RTSuccessful(i) = firstLickTime - secondToneTime;
            else %FOR TEST PROTOCOL 
                % If the first lick time is NOT after the second tone, store it in the unsuccessful RT array
                RTUnsuccessful(i) = firstLickTime - firstToneTime;
            end  
        else %FOR TRAIN PROTOCOL
            % If the first lick time is NOT after the second tone, store it in the unsuccessful RT array
            RTUnsuccessful(i) = firstLickTime - firstToneTime;
        end            
    else %FOR TRAIN PROTOCOL
        FirstTone(i) = NaN;
    end
end

%Remove the First Trial from PSTH Matrix
psthMatrix(1, :, :) = [];

% Define unsuccessful trials based on RT
UnsuccessfulTrialLogical = isnan(RTSuccessful) | RTSuccessful > RTThresholdmax | RTSuccessful < RTThresholdmin;
UnsuccessfulTrialLogical2 = isnan(RTUnsuccessful) | RTUnsuccessful < 0;

% Use the opposite logic for successful trials (i.e., NOT unsuccessful)
SuccessfulTrialLogical = ~UnsuccessfulTrialLogical;
SuccessfulTrialLogical2 = ~UnsuccessfulTrialLogical2;

SuccessfulTrial = TrialStartAdj(SuccessfulTrialLogical);
FirstToneSuccessful = FirstTone(SuccessfulTrialLogical);
FirstToneTimesSuccess = SuccessfulTrial + FirstToneSuccessful;  % Time relative to photometry start
RTToneSuccessful = RTSuccessful(SuccessfulTrialLogical);
RTToneUnsuccessful = RTUnsuccessful(SuccessfulTrialLogical2);
ITISuccessful = ITI(SuccessfulTrialLogical);
ITIUnsuccessful = ITI(SuccessfulTrialLogical2);
FirstToneTimesWithoutSecond = TrialStartAdj(SuccessfulTrialLogical2) + FirstTone(SuccessfulTrialLogical2);

% Extract successful PSTH data
for a = 1:numfibers
    psthSuccessArray = NaN(numel(FirstToneTimesSuccess), (tPrev+tPost+1));
    for x = 1:numel(FirstToneTimesSuccess)
        [~, thisTime] = min(abs(time-FirstToneTimesSuccess(x)));
        
        % Check for condition where (thisTime+tPost) exceeds array bounds
        if (x == numel(FirstToneTimesSuccess)) && ((thisTime+tPost) > length(sig_normZ))
            FirstToneTimesSuccess(end) = [];
            RTToneSuccessful(end) = [];
            RTSuccessful(end) = [];
            ITISuccessful(end) = [];
            warning('The last trial has been removed due to exceeding array bounds.');
            psthSuccessArray(end, :) = []; % Remove the last row of the current psthSuccessArray
            continue; % Skip this iteration
        end
        
        psthSuccessArray(x, :) = sig_normZ(a, thisTime-tPrev:thisTime+tPost); %normalize to first tone
    end
    psthSuccessMatrix(:, :, a) = psthSuccessArray;
end

% Generate the psth matrix for trials that have the first tone but not the second tone
for a = 1:numfibers
    psthWithoutSecondArray = NaN(numel(FirstToneTimesWithoutSecond), (tPrev+tPost+1));  
    for x = 1:numel(FirstToneTimesWithoutSecond)
        [~, thisTime] = min(abs(time - FirstToneTimesWithoutSecond(x)));
        
        % Check for condition where (thisTime+tPost) exceeds array bounds
        if (x == numel(FirstToneTimesWithoutSecond)) && ((thisTime+tPost) > length(sig_normZ))
            FirstToneTimesWithoutSecond(end) = [];
            RTToneUnsuccessful(end) = [];
            RTUnsuccessful(end) = [];
            ITIUnsuccessful(end) = [];
            warning('The last trial has been removed due to exceeding array bounds.');
            psthWithoutSecondArray(end, :) = []; % Remove the last row of the current psthWithoutSecondArray
            continue; % Skip this iteration
        end
        
        psthWithoutSecondArray(x, :) = sig_normZ(a, thisTime-tPrev:thisTime+tPost); %normalize to first tone
    end
    psthWithoutSecondMatrix(:, :, a) = psthWithoutSecondArray;
end

end
