% CCDT_feature_selection
% robust feature selection with bootstrapped regression using random 80% of trials
% MH 09/2022
clc; clear; close all
warning off

% parameters
db = CCDTdatabase;
stsubj = 1; % start subject
sigPThresh = 0.05; %pvalue threshold
nReps = 1000; % # repetitions for bootstrap
trainPerc = .8; % fraction of trials to train on
fbands = [3 12; 12 30; 30 55; 70 110]; 
Nsubj = size(db,1); Nbands = size(fbands,1)+1; Nfreq = height(fbands);

% feature data info
pdir = 'test_output\'; % feature directory
pfname = 'power_features.mat'; % power feature file name 
gfname = 'graph_features.mat'; % graph communicability feature file name

% output file info
savedir = 'test_output\'; %output directory
savefname = 'feature_selection'; %output file name

% load NFstruct for feature data and behavStruct for RT info
load([pdir pfname],"NFstruct"); 
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","behavStruct");
NFstructQ = NFstruct;
clear NFstruct;

% initiate empty structures
trainIdx = cell(Nsubj,1);
testIdx = cell(Nsubj,1);
sigPow = cell(Nsubj,3);
sigCom = cell(Nsubj,3);
LTdat = cell(Nsubj,2); %size = (Nsubj,2,4) 2nd dimension: power features, qexp features, 3rd dimension index: 1 = slope estimate, 2 = ci 1, 3 = ci 2, 4 = pval

for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    
    % organize pow/com and RT data for isubj
    for iFreq = 1:Nfreq
        comDat(:,:,iFreq) = NFstructQ(isubj).fbands(iFreq).NFqexp_nodal;
        powDat(:,:,iFreq) =  NFstructP(isubj).fbands(iFreq).NFpow;

        NchQ = size(comDat,1);
        NchP = size(powDat,1);
    end
    
    % remove trials with bad RTs from analysis
    cRT = behavStruct(isubj).vRT;
    goodTrl = ~(cRT<50 | cRT>999); % RT<50 too fast to be physiologically realistic response to cue
    nTrl = length(find(goodTrl));
    trainTrls = floor(nTrl*trainPerc);
    cRT = cRT(goodTrl);
    comDat = comDat(:,goodTrl,:);
    powDat = powDat(:,goodTrl,:);

    %bootstrap loop
    trainIdx_temp = zeros(trainTrls,nReps);
    testIdx_temp = zeros(nTrl-trainTrls,nReps);
    powPval = zeros(NchP, length(fbands), nReps);
    comPval = zeros(NchQ, length(fbands), nReps);
    powBval = zeros(NchP, length(fbands), nReps);
    comBval = zeros(NchQ, length(fbands), nReps);
    parfor i = 1:nReps
        if ~mod(i,100)
            disp(i) %print every 100 reps
        end

        % choose random 80% of trials
        idxShuffle = randperm(nTrl); %random trial order
        train = idxShuffle(1:trainTrls);
        trainIdx_temp(:,i) = train;
        testIdx_temp(:,i) = idxShuffle(trainTrls+1:end);

        % allocate local results
        powBval_temp = zeros(NchP, length(fbands));
        powPval_temp = zeros(NchP, length(fbands));
        comBval_temp = zeros(NchQ, length(fbands));
        comPval_temp = zeros(NchQ, length(fbands));
        
        % univariate feature selection
        for iChan = 1:NchP % power features
            for iFreq = 1:length(fbands)
                [b, stats] = robustfit(powDat(iChan,train,iFreq),cRT(train));
                powPval_temp(iChan,iFreq) = stats.p(2);
                powBval_temp(iChan,iFreq) = b(2);
%                 [bpow(iChan,iFreq,:),statspow(iChan,iFreq,:)] = robustfit(powDat(iChan,train,iFreq),cRT(train));
%                 powPval(iChan,iFreq,i) = statspow(iChan,iFreq).p(2);
                % LTdat 3rd dimension index: 1 = slope estimate, 2 = ci 1, 3 = ci 2, 4 = pval
            end
        end
        for iChan = 1:NchQ % communicabililty features
            for iFreq = 1:length(fbands)
                [b, stats] = robustfit(comDat(iChan,train,iFreq),cRT(train));
                comPval_temp(iChan,iFreq) = stats.p(2);
                comBval_temp(iChan,iFreq) = b(2);

%                 [bcom(iChan,iFreq,:),statscom(iChan,iFreq,:)] = robustfit(comDat(iChan,train,iFreq),cRT(train));
%                 LTdat{isubj,2}(iChan,iFreq,:) = [bcom(iChan,iFreq,2) bcom(iChan,iFreq,2)-1.96*statscom(iChan,iFreq).se(2) bcom(iChan,iFreq,2)+1.96*statscom(iChan,iFreq).se(2) statscom(iChan,iFreq).p(2)];
%                 comPval(iChan,iFreq,i) = statscom(iChan,iFreq).p(2);
                % LTdat 3rd dimension index: 1 = slope estimate, 2 = ci 1, 3 = ci 2, 4 = pval
            end
        end
        powBval(:,:,i) = powBval_temp;
        powPval(:,:,i) = powPval_temp;
        comBval(:,:,i) = comBval_temp;
        comPval(:,:,i) = comPval_temp;
    end

    %sigPow{:,1} = binary significance. sigPow{:,2} = pvalues. sigPow{:,3} = slope values
    sigPow{isubj,1} = powPval <= sigPThresh;
    sigCom{isubj,1} = comPval <= sigPThresh;
    sigPow{isubj,2} = powPval;
    sigCom{isubj,2} = comPval;
    sigPow{isubj,3} = powBval;
    sigCom{isubj,3} = comBval;
    trainIdx{isubj} = trainIdx_temp;
    testIdx{isubj} = testIdx_temp;
    
    clear comBval* powBval* goodTrl bcom* bpow* statscom* statspow* comDat* powDat* comPval* powPval* trainIdx_temp testIdx_temp
end

% save feature selection info from bootstrapped regressions
% used below and in CCDT_Classifier
save([savedir savefname '.mat'],"trainIdx","testIdx","sigPow","sigCom")
