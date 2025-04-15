% updated to use with data processed with Vivek's processing scripts, not Brandon's
clc; clear; close all
warning off

% parameters
db = CCDTdatabase;
stsubj = 1; % start subject
sigPThresh = 0.05; %pvalue threshold
newPermutation = 0; % do univariate feature selection with new permutations of 80% of trials
nReps = 1000; % # repetitions for bootstrap (if newPermutation = 1)
trainPerc = .8;
fbands = [3 12; 12 30; 30 55; 70 110]; 
Nsubj = size(db,1); Nbands = size(fbands,1)+1; Nfreq = height(fbands);

% feature data info
pdir = ''; % feature directory
pfname = ''; % power feature file name
gfname = ''; % graph communicability feature file name

% behavioral data info
rtdir = '/mnt/sdb1/CCDT/orig_procData/'; % behavioral data directory
rtfnm = 'graphRT_sIall_500mspre'; %use graph metric versus PLV metric defined here

savedir = ''; %output directory
savefname = ''; %output file name
loadfname = ''; % file to load indices from if newPermutation = 0

%load NFstruct data and pvals from feature selection with 100% of trials
load([pdir pfname],"NFstruct","LTpow"); %
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","LTqexp");
NFstructQ = NFstruct;
clear NFstruct;
load([rtdir rtfname],"behavStruct"); %load behavStruct RT data

if ~newPermutation
    %loads testIdx, trainIdx
    load([savedir loadfname '.mat'],"trainIdx","testIdx")  
end
if newPermutation
    trainIdx = cell(Nsubj,1);
    testIdx = cell(Nsubj,1);
end

nullShuffleRT = cell(Nsubj,1);
nullRandData = cell(Nsubj,4);

for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    
    % organize pow/com and RT data for isubj
    for iFreq = 1:Nfreq
        comDat(:,:,iFreq) = NFstructQ(isubj).fbands(iFreq).NFqexp_nodal;
        powDat(:,:,iFreq) =  NFstructP(isubj).fbands(iFreq).NFpow;

        NchQ = size(comDat,1);
        NchP = size(powDat,1);
    end
    
    cRT = behavStruct(isubj).vRT;
    goodTrl = ~(cRT<50 | cRT>999);
    nTrl = length(find(goodTrl));
    trainTrls = floor(nTrl*trainPerc);

    cRT = cRT(goodTrl);
    comDat = comDat(:,goodTrl,:);
    powDat = powDat(:,goodTrl,:);
  
    % do analysis

    % set up null comparisons:
    % random permutations of RT 
    nullRT = cRT(randperm(length(cRT)));
    nullShuffleRT{isubj} = nullRT;
    % add random noise to data
    nullPowDat = powDat.*randn(size(powDat));
    nullRandData{isubj,1} = nullPowDat;
    nullComDat = comDat.*randn(size(comDat));
    nullRandData{isubj,2} = nullComDat;
    nullPowDatAdd = powDat+randn(size(powDat));
    nullRandData{isubj,3} = nullPowDatAdd;
    nullComDatAdd = comDat+randn(size(comDat));
    nullRandData{isubj,4} = nullComDatAdd;
    if newPermutation
        trainIdx{isubj} = zeros(trainTrls);
        testIdx{isubj} = zeros(nTrl-trainTrls);
    end
    %bootstrap loop
    for i = 1:nReps
        if ~mod(i,100)
            disp(i) %print every 100 reps
        end

        if newPermutation
            % choose random 80% of trials
            idxShuffle = randperm(nTrl); %random trial order
            train = idxShuffle(1:trainTrls);
            trainIdx{isubj}(i,:) = train;
            testIdx{isubj}(i,:) = idxShuffle(trainTrls+1:end);
        else
            train = trainIdx{isubj}(i,:);
        end
        % univariate feature selection
        for iChan = 1:NchP
            for iFreq = 1:length(fbands)
                % behavioral data
                [bpow(iChan,iFreq,:),statspow(iChan,iFreq,:)] = robustfit(powDat(iChan,train,iFreq),cRT(train));
                LTdat{isubj,1}(iChan,iFreq,:) = [bpow(iChan,iFreq,2) bpow(iChan,iFreq,2)-1.96*statspow(iChan,iFreq).se(2) bpow(iChan,iFreq,2)+1.96*statspow(iChan,iFreq).se(2) statspow(iChan,iFreq).p(2)];
                powPval(iChan,iFreq,i) = statspow(iChan,iFreq).p(2);
                % LTdat 3rd dimension index: 1 = estimate, 2 = ci 1, 3 = ci 2, 4 = pval

%               % shuffled RT null
                [bpowNRT(iChan,iFreq,:),statspowNRT(iChan,iFreq,:)] = robustfit(powDat(iChan,train,iFreq),nullRT(train));
                LTdatNRT{isubj,1}(iChan,iFreq,:) = [bpowNRT(iChan,iFreq,2) bpowNRT(iChan,iFreq,2)-1.96*statspowNRT(iChan,iFreq).se(2) bpowNRT(iChan,iFreq,2)+1.96*statspowNRT(iChan,iFreq).se(2) statspowNRT(iChan,iFreq).p(2)];
                powPvalNRT(iChan,iFreq,i) = statspowNRT(iChan,iFreq).p(2);
                
            end
        end
        for iChan = 1:NchQ
            for iFreq = 1:length(fbands)
                % behavioral data
                [bcom(iChan,iFreq,:),statscom(iChan,iFreq,:)] = robustfit(comDat(iChan,train,iFreq),cRT(train));
                LTdat{isubj,2}(iChan,iFreq,:) = [bcom(iChan,iFreq,2) bcom(iChan,iFreq,2)-1.96*statscom(iChan,iFreq).se(2) bcom(iChan,iFreq,2)+1.96*statscom(iChan,iFreq).se(2) statscom(iChan,iFreq).p(2)];
                comPval(iChan,iFreq,i) = statscom(iChan,iFreq).p(2);
                % LTdat 3rd dimension index: 1 = estimate, 2 = ci 1, 3 = ci 2, 4 = pval

                % shuffled RT null
                [bcomNRT(iChan,iFreq,:),statscomNRT(iChan,iFreq,:)] = robustfit(comDat(iChan,train,iFreq),nullRT(train));
                LTdatNRT{isubj,2}(iChan,iFreq,:) = [bcomNRT(iChan,iFreq,2) bcomNRT(iChan,iFreq,2)-1.96*statscomNRT(iChan,iFreq).se(2) bcomNRT(iChan,iFreq,2)+1.96*statscomNRT(iChan,iFreq).se(2) statscomNRT(iChan,iFreq).p(2)];
                comPvalNRT(iChan,iFreq,i) = statscomNRT(iChan,iFreq).p(2);
            end
        end
        powBval(:,:,i) = bpow(:,:,2);
        comBval(:,:,i) = bcom(:,:,2);
        powBvalNRT(:,:,i) = bpowNRT(:,:,2);
        comBvalNRT(:,:,i) = bcomNRT(:,:,2);
    end

    %sigPow{:,1} = binary significance. sigPow{:,2} = pvalues. sigPow{:,3} = slope values
    sigPow{isubj,1} = powPval <= sigPThresh;
    sigCom{isubj,1} = comPval <= sigPThresh;
    sigPow{isubj,2} = powPval;
    sigCom{isubj,2} = comPval;
    sigPow{isubj,3} = powBval;
    sigCom{isubj,3} = comBval;

    sigPowNRT{isubj,1} = powPvalNRT <= sigPThresh;
    sigComNRT{isubj,1} = comPvalNRT <= sigPThresh;
    sigPowNRT{isubj,2} = powPvalNRT;
    sigComNRT{isubj,2} = comPvalNRT;
    sigPowNRT{isubj,3} = powBvalNRT;
    sigComNRT{isubj,3} = comBvalNRT;
    
    clear comBval* powBval* goodTrl bcom* bpow* statscom* statspow* comDat* powDat* comPval* powPval*
end

if ~(exist('testIdx','var') == 1)
    load([savedir loadfname '.mat'],"trainIdx","testIdx")  
end
save([savedir savefname '.mat'],"trainIdx","testIdx","nullShuffleRT","sigPow","sigCom","sigPowNRT","sigComNRT")
