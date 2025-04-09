% updated to use with data processed with Vivek's processing scripts, not Brandon's
clc; clear; close all
warning off
db = CCDTdatabase;

% parameters
stsubj = 1; % start subject
rrfNames = {'bipolar' 'commonAvg' 'commonMed'};
fbandNames = ["T/A" "B" "loG" "hiG"];
pdir = '/mnt/sdb1/CCDT/procData072023/';
% pfname = 'powRT_071023_Bipolar_BR_zAll';
% pfname = 'powRT_082823_CAR_BR_zAll';
% pfname = 'powRT_090923_BP_intra_MH_zAll';
% pfname = 'powRT_090923_CAR_intra_MH_zAll';
% pfname = 'powRT_101123_Bipolar_zAll'; %main data to use
% pfname = 'powRT_103123_Bipolar_zAll_null2';
% pfname = 'powRT_031324_CAR_hg70-110'; % high gamma range 70-100Hz
% pfname = 'powRT_090924_CAR_hg70110';
pfname = 'powRT_040425_CAR_hg70110';
% pfname = 'powRT_090924_CAR_intra_hg70110';
% pfname = 'powRT_043024_BP_remove60HZoutl'; % removed outliers w/ high 60hz power after processing
% pfname = 'powRT_050624_BP_intra_hg70-110_removeOutl';
% gfname = 'graphRT_090923_CAR_intra_MH';
% gfname = 'graphRTpli_083123_CAR_BR';
% gfname = 'graphRT_071823_CAR_BR'; %main data to use
% gfname = 'graphRT_103123_CAR_null2'; %random circular shift
% gfname = 'graphRT_020424_CAR_PLV'; % PLV instead of qexp as graph metric
% gfname = 'graphRT_031024_CAR_hg70-110'; 
gfname = 'graphRT_040425_CAR_hg70110'; %
% gfname = 'graphRT_090324_CAR_hg70110'; %pretrial qexp
% gfname = 'graphRT_031025_CAR_PLV_hg70110'; %pretrial plv
% gfname = 'graphRT_031025_CAR_PLVSS1_hg70110'; %pretrial PLV spatial smoothing
% gfname = 'graphRT_090324_CAR_intra_hg70110'; %intratrial qexp
% gfname = 'graphRT_050624_CAR_intra_hg70-110'; % high gamma range 70-100Hz
% gfname = 'graphRT_102423_CAR_intra';
rtdir = '/mnt/sdb1/CCDT/CCDT Data/CCDT/procData/';
rtfname = 'graphRT_sIall_500mspre';
sigPThresh = 0.05; %pvalue threshold
doZ = 0; % zscore data?
zType = [1 2]; % 1 = trial, 2 = electrode, 3 = freq, 'all' = whole matrix
newPermutation = 0; % do univariate feature selection with new permutations of 80% of trials
% percThresh = .5; % threshold of what % of repetitions to consider a feature selected (if newPermutation = 0)
nReps = 1000; % # repetitions for bootstrap (if newPermutation = 1)
trainPerc = .8;
savedir = '/mnt/sdb1/CCDT/CCDTscripts/80-20 dat/';
savefname = '80-20-dat-pCAR-qCAR-hg70110-040425'; % file to save 
loadfname = '80-20-dat-pBP-qCAR-noZ'; % file to load indices from if newPermutation = 0
% loadfname = '80-20-dat-pBP-qCAR-intra'; %train/test idx for intratrial data
fbands = [3 12; 12 30; 30 55; 70 110]; %[3 12; 12 30; 30 55; 70 150];
Nsubj = size(db,1); Nbands = size(fbands,1)+1; Nfreq = height(fbands);

%load NFstruct data and pvals from feature selection with 100% of trials
load([pdir pfname],"NFstruct","LTpow"); %
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","LTqexp");
NFstructQ = NFstruct;
clear NFstruct;
%load behavStruct RT data
load([rtdir rtfname],"behavStruct");

if ~newPermutation
    %loads sigCom, sigPow, sigStr, testIdx, trainIdx
    load([savedir loadfname '.mat'],"trainIdx","testIdx")  
end

nullShuffleRT = cell(Nsubj,1);
nullRandData = cell(Nsubj,4); %{nullPowDat, nullComDat, nullPowDatAdd, nullComDatAdd}
if newPermutation
    trainIdx = cell(Nsubj,1);
    testIdx = cell(Nsubj,1);
end

for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    
    % organize pow/com and RT data for isubj
    for iFreq = 1:Nfreq
        comDat(:,:,iFreq) = NFstructQ(isubj).fbands(iFreq).NFqexp_nodal;
        powDat(:,:,iFreq) =  NFstructP(isubj).fbands(iFreq).NFpow;

        NchQ = size(comDat,1);
        NchP = size(powDat,1);
    end
    if doZ
        powDat = zscore(powDat,[],zType);
        comDat = zscore(comDat,[],zType);
    end

    cRT = behavStruct(isubj).vRT;
    goodTrl = ~(cRT<50 | cRT>999);
    nTrl = length(find(goodTrl));
    trainTrls = floor(nTrl*trainPerc);

    cRT = cRT(goodTrl);
    comDat = comDat(:,goodTrl,:);
    powDat = powDat(:,goodTrl,:);
  
    % do analysis
%     if newPermutation
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

%                     % shuffled RT null
                    [bpowNRT(iChan,iFreq,:),statspowNRT(iChan,iFreq,:)] = robustfit(powDat(iChan,train,iFreq),nullRT(train));
                    LTdatNRT{isubj,1}(iChan,iFreq,:) = [bpowNRT(iChan,iFreq,2) bpowNRT(iChan,iFreq,2)-1.96*statspowNRT(iChan,iFreq).se(2) bpowNRT(iChan,iFreq,2)+1.96*statspowNRT(iChan,iFreq).se(2) statspowNRT(iChan,iFreq).p(2)];
                    powPvalNRT(iChan,iFreq,i) = statspowNRT(iChan,iFreq).p(2);
                    
%                     % random noise in data null (multiply)
%                     [bpowNdat(iChan,iFreq,:),statspowNdat(iChan,iFreq,:)] = robustfit(nullPowDat(iChan,train,iFreq),cRT(train));
%                     LTdatNdat{isubj,1}(iChan,iFreq,:) = [bpowNdat(iChan,iFreq,2) bpowNdat(iChan,iFreq,2)-1.96*statspowNdat(iChan,iFreq).se(2) bpowNdat(iChan,iFreq,2)+1.96*statspowNdat(iChan,iFreq).se(2) statspowNdat(iChan,iFreq).p(2)];
%                     powPvalNdat(iChan,iFreq,i) = statspowNdat(iChan,iFreq).p(2);
% 
%                      % random noise in data null (add)
%                     [bpowNdatAdd(iChan,iFreq,:),statspowNdatAdd(iChan,iFreq,:)] = robustfit(nullPowDatAdd(iChan,train,iFreq),cRT(train));
%                     LTdatNdatAdd{isubj,1}(iChan,iFreq,:) = [bpowNdatAdd(iChan,iFreq,2) bpowNdatAdd(iChan,iFreq,2)-1.96*statspowNdatAdd(iChan,iFreq).se(2) bpowNdatAdd(iChan,iFreq,2)+1.96*statspowNdatAdd(iChan,iFreq).se(2) statspowNdatAdd(iChan,iFreq).p(2)];
%                     powPvalNdatAdd(iChan,iFreq,i) = statspowNdatAdd(iChan,iFreq).p(2);
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

%                     % random noise in data null (multiply)
%                     [bcomNdat(iChan,iFreq,:),statscomNdat(iChan,iFreq,:)] = robustfit(nullComDat(iChan,train,iFreq),cRT(train));
%                     LTdatNdat{isubj,2}(iChan,iFreq,:) = [bcomNdat(iChan,iFreq,2) bcomNdat(iChan,iFreq,2)-1.96*statscomNdat(iChan,iFreq).se(2) bcomNdat(iChan,iFreq,2)+1.96*statscomNdat(iChan,iFreq).se(2) statscomNdat(iChan,iFreq).p(2)];
%                     comPvalNdat(iChan,iFreq,i) = statscomNdat(iChan,iFreq).p(2);

%                     % random noise in data null (multiply)
%                     [bcomNdatAdd(iChan,iFreq,:),statscomNdatAdd(iChan,iFreq,:)] = robustfit(nullComDatAdd(iChan,train,iFreq),cRT(train));
%                     LTdatNdatAdd{isubj,2}(iChan,iFreq,:) = [bcomNdatAdd(iChan,iFreq,2) bcomNdatAdd(iChan,iFreq,2)-1.96*statscomNdatAdd(iChan,iFreq).se(2) bcomNdatAdd(iChan,iFreq,2)+1.96*statscomNdatAdd(iChan,iFreq).se(2) statscomNdatAdd(iChan,iFreq).p(2)];
%                     comPvalNdatAdd(iChan,iFreq,i) = statscomNdatAdd(iChan,iFreq).p(2);
                end
            end
            powBval(:,:,i) = bpow(:,:,2);
            comBval(:,:,i) = bcom(:,:,2);
            powBvalNRT(:,:,i) = bpowNRT(:,:,2);
            comBvalNRT(:,:,i) = bcomNRT(:,:,2);
%             powBvalNdat(:,:,i) = bpowNdat(:,:,2);
%             comBvalNdat(:,:,i) = bcomNdat(:,:,2);
%             powBvalNdatAdd(:,:,i) = bpowNdatAdd(:,:,2);
%             comBvalNdatAdd(:,:,i) = bcomNdatAdd(:,:,2);
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
% 
%         sigPowNdat{isubj,1} = powPvalNdat <= sigPThresh;
%         sigComNdat{isubj,1} = comPvalNdat <= sigPThresh;
%         sigPowNdat{isubj,2} = powPvalNdat;
%         sigComNdat{isubj,2} = comPvalNdat;
%         sigPowNdat{isubj,3} = powBvalNdat;
%         sigComNdat{isubj,3} = comBvalNdat;
        
%         sigPowNdatAdd{isubj,1} = powPvalNdatAdd <= sigPThresh;
%         sigComNdatAdd{isubj,1} = comPvalNdatAdd <= sigPThresh;
%         sigPowNdatAdd{isubj,2} = powPvalNdatAdd;
%         sigComNdatAdd{isubj,2} = comPvalNdatAdd;
%         sigPowNdatAdd{isubj,3} = powBvalNdatAdd;
%         sigComNdatAdd{isubj,3} = comBvalNdatAdd;
        
        clear comBval* powBval* 
        
%     else
%         % use permutations loaded from 80-20-dat
%         % feature selection already done
%         % comparing selected features using 100% of data vs 80% of data
%         
%         sigPowPerc = sum(sigPow{isubj,1},3)/nReps;
%         sigComPerc = sum(sigCom{isubj,1},3)/nReps;
% 
%         [sigElecPow,sigFreqPow] = find(sigPowPerc>percThresh);
%         chlabsPow = NFstructP(isubj).gchlbl(sigElecPow);
% %             chlabsPow = elecDat.label(sigElecPow); %electrode numbers
%         flabsPow = fbandNames(sigFreqPow);
%         
%         [sigElecCom,sigFreqCom] = find(sigComPerc>percThresh);
%         chlabsCom = NFstructQ(isubj).gchlbl(sigElecCom);
% %             chlabsCom = elecDat.label(sigElecCom); %electrode numbers
%         flabsCom = fbandNames(sigFreqCom);
% 
%         % create 80-20 table
%         sigContacts8020{isubj} = table;
%         sigContacts8020{isubj}.Metric = categorical([repmat("pow", [length(sigElecPow),1]); repmat("com", [length(sigElecCom),1])]);
%         sigContacts8020{isubj}.Channel = categorical([chlabsPow;chlabsCom]);
%         sigContacts8020{isubj}.Fband = categorical([flabsPow,flabsCom]');
% %             sigContacts8020{isubj}.aBNlab = categorical([elecDat.aBNlab(sigElecPow); elecDat.aBNlab(sigElecCom)]);
% %             sigContacts8020{isubj}.aEVElab = categorical([elecDat.aEVElab(sigElecPow); elecDat.aEVElab(sigElecCom)]);
% %             sigContacts8020{isubj}.BNeveG = categorical([elecDat.BNeveG(sigElecPow); elecDat.BNeveG(sigElecCom)]);
% 
%         ch = [chlabsPow;chlabsCom];
%         f = [flabsPow,flabsCom]';
%         m = [repmat("pow", [length(sigElecPow),1]); repmat("com", [length(sigElecCom),1])];
% %             groupLab = string([elecDat.BNeveG(sigElecPow); elecDat.BNeveG(sigElecCom)]);
%         undsc = repmat("-", [length(ch),1]);
%         sigContacts8020{isubj}.NodeFeature = categorical(strcat(ch,undsc,m,undsc,f));
% %             sigContacts8020{isubj}.AnatNodeFeature = categorical(strcat(groupLab,undsc,m,undsc,f));
%         
% 
%         sigNodesPow = LTpow{isubj}(:,:,4)<=sigPThresh;
%         [sigElecPow100,sigFreqPow100] = find(sigNodesPow);
%         chlabsPow100 = NFstructP(isubj).gchlbl(sigElecPow100); %electrode numbers
%         flabsPow100 = fbandNames(sigFreqPow100);
%         
%         sigNodesCom = LTqexp{isubj}(:,:,4)<=sigPThresh;
%         [sigElecCom100,sigFreqCom100] = find(sigNodesCom);
%         chlabsCom100 = NFstructQ(isubj).gchlbl(sigElecCom100); %electrode numbers
%         flabsCom100 = fbandNames(sigFreqCom100);
% 
%         % create 100 table
%         sigContacts100{isubj} = table;
%         sigContacts100{isubj}.Metric = categorical([repmat("pow", [length(sigElecPow100),1]); repmat("com", [length(sigElecCom100),1])]);
%         sigContacts100{isubj}.Channel = categorical([chlabsPow100;chlabsCom100]);
%         sigContacts100{isubj}.Fband = categorical([flabsPow100,flabsCom100]');
% %             sigContacts100{isubj}.aBNlab = categorical([elecDat.aBNlab(sigElecPow100); elecDat.aBNlab(sigElecCom100)]);
% %             sigContacts100{isubj}.aEVElab = categorical([elecDat.aEVElab(sigElecPow100); elecDat.aEVElab(sigElecCom100)]);
% %             sigContacts100{isubj}.BNeveG = categorical([elecDat.BNeveG(sigElecPow100); elecDat.BNeveG(sigElecCom100)]);
% 
%         ch = [chlabsPow100;chlabsCom100];
%         f = [flabsPow100,flabsCom100]';
%         m = [repmat("pow", [length(sigElecPow100),1]); repmat("com", [length(sigElecCom100),1])];
% %             groupLab = string([elecDat.BNeveG(sigElecPow100); elecDat.BNeveG(sigElecCom100)]);
%         undsc = repmat("-", [length(ch),1]);
%         sigContacts100{isubj}.NodeFeature = categorical(strcat(ch,undsc,m,undsc,f));
%             sigContacts100{isubj}.AnatNodeFeature = categorical(strcat(groupLab,undsc,m,undsc,f));
%         
%     end

    clear goodTrl bcom* bpow* statscom* statspow* comDat* powDat* comPval* powPval*
        
end
%%
% fix sigCom/sigPow dimensions
% for i = 1:Nsubj
%     %get # of good channels
%     NchQ = size(NFstructQ(i).gch,1);
%     NchP = size(NFstructP(i).gch,1);
% 
%     for j = 1:3
%         sigCom{i,j} = sigCom{i,j}(1:NchQ,:,:);
%         sigComNRT{i,j} = sigComNRT{i,j}(1:NchQ,:,:);
%         sigComNdat{i,j} = sigComNdat{i,j}(1:NchQ,:,:);
%         sigPow{i,j} = sigPow{i,j}(1:NchP,:,:);
%         sigPowNRT{i,j} = sigPowNRT{i,j}(1:NchP,:,:);
%         sigPowNdat{i,j} = sigPowNdat{i,j}(1:NchP,:,:);
%     end
% end
%%
if ~(exist('testIdx','var') == 1)
    load([savedir loadfname '.mat'],"trainIdx","testIdx")  
end
% % save data from univariate feature selection
% save([savedir savefname '.mat'],"trainIdx","testIdx","nullShuffleRT","nullRandData","sigPow","sigCom","sigPowNRT","sigComNRT","sigPowNdat","sigComNdat","sigPowNdatAdd","sigComNdatAdd")
% save([savedir savefname '.mat'],"trainIdx","testIdx","nullRandData","sigPowNdatAdd","sigComNdatAdd")
save([savedir savefname '.mat'],"trainIdx","testIdx","nullShuffleRT","sigPow","sigCom","sigPowNRT","sigComNRT")
% sigComNdat = sigCom;
% sigPowNdat = sigPow;
% save([savedir savefname '.mat'],"trainIdx","testIdx","sigPowNdat","sigComNdat")

% save([savedir savefname '.mat'],"trainIdx","testIdx","sigPow","sigCom") %normal option

