% function CCDTnodeAnalyzeNew_MH
%   MH 08/2023
% use for generating fPv/sPv/fQv/sQv etc 

clear; close all; clc

db = CCDTdatabase;
Nsubj = height(db);

saveFile = 1;
savedir = ''; %output directory
savefname = ''; %output file name
percThresh = .3; % threshold w/ highest mean AUC from SVM analysis

%load behavioral data
rtdir = '/mnt/sdb1/CCDT/orig_procData/'; % behavioral data directory
rtfnm = 'graphRT_sIall_500mspre'; %use graph metric versus PLV metric defined here
load([rtdir rtfnm '.mat'],'behavStruct'); %get RT and fast/slow trial info from this

%load features
pdir = ''; % feature directory
pfname = ''; % power feature file name
gfname = ''; % graph communicability feature file name
load([pdir pfname],"NFstruct","LTpow"); 
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","LTqexp");
NFstructQ = NFstruct;
clear NFstruct;

% load selected features
sdir = ''; %feature selection directory
fname8020 = ''; %feature selection file name
load([sdir fname8020],'sigPow','sigCom')

% get features above selection threshold
for isubj = 1:Nsubj
    sigNodesPow{isubj,1} = (sum(sigPow{isubj,1},3)/1000 > percThresh);
    sigNodesCom{isubj,1} = (sum(sigCom{isubj,1},3)/1000 > percThresh);
end


%% analysis section
fPv = []; sPv = []; %fast and slow power features
fIDp = []; ssIDp = []; %frequency band of power features
fQv = []; sQv = []; %fast and slow communicability features
aPv = []; aQv = []; %all power and communicability features
fIDq = []; ssIDq = []; %frequency band of communicability features
gchP = []; gchQ = []; %channel labels of power and communicability features
for isubj = 1:Nsubj
    disp(['Subject: ' num2str(isubj)])

    % get iF, iS - fast and slow trial indices
    RT = behavStruct(isubj).vRT;
    iF = behavStruct(isubj).ifast;
    iS = behavStruct(isubj).islow;

    for j = 1:4 %fbands, reformat features
        nComDat(:,:,j) = NFstructQ(isubj).fbands(j).NFqexp_nodal;
        nPowDat(:,:,j) = NFstructP(isubj).fbands(j).NFpow;
    end

    sigFeatP = sigNodesPow{isubj}; %nElec x nFreq
    nElecP = size(nPowDat,1);
    nTrl = size(nPowDat,2);
    nFreq = size(nPowDat,3);
    sigFeatQ = sigNodesCom{isubj};
    nElecCAR = size(nComDat,1);
    
    %get selected features
    for i = 1:nFreq
        powFeaturesZ = zscore(nPowDat(:,:,i),1,'all'); % zscore over all channels/trials
        sigPowFeaturesZ = powFeaturesZ(sigFeatP(:,i),:); % get just selected channels
        fPv = vertcat(fPv,mean(sigPowFeaturesZ(:,iF),2)); %average over fast trials
        sPv = vertcat(sPv,mean(sigPowFeaturesZ(:,iS),2)); %average over slow trials
        aPv = vertcat(aPv,mean(sigPowFeaturesZ(:,:),2)); %average over all trials
        if ~isempty(sigPowFeaturesZ)
            fIDp = vertcat(fIDp,i*ones(height(sigPowFeaturesZ),1));
            ssIDp = vertcat(ssIDp,isubj*ones(height(sigPowFeaturesZ),1));
            gchP = vertcat(gchP,NFstructP(isubj).gchlbl(sigFeatP(:,i)));
        end
        comFeaturesZ = zscore(nComDat(:,:,i),1,'all');
        sigComFeaturesZ = comFeaturesZ(sigFeatQ(:,i),:);
        fQv = vertcat(fQv,mean(sigComFeaturesZ(:,iF),2)); %average over fast trials
        sQv = vertcat(sQv,mean(sigComFeaturesZ(:,iS),2)); %average over slow trials
        aQv = vertcat(aQv,mean(sigComFeaturesZ(:,:),2)); %average over all trials
        if ~isempty(sigComFeaturesZ)
            fIDq = vertcat(fIDq,i*ones(height(sigComFeaturesZ),1));
            ssIDq = vertcat(ssIDq,isubj*ones(height(sigComFeaturesZ),1));
            gchQ = vertcat(gchQ,NFstructQ(isubj).gchlbl(sigFeatQ(:,i)));
        end
    end
    clear nComDat nPowDat
end

if saveFile 
    save([savedir savefname],'fPv','sPv','fIDp','ssIDp','gchP','fQv','sQv','aQv','aPv','fIDq','ssIDq','gchQ')
end