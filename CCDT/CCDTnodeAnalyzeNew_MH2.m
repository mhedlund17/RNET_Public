% function CCDTnodeAnalyzeNew_MH
%   MH 08/2023
% use for generating fPv/sPv/fQv/sQv etc for data we processed like
% original data with features selected from 80% of trials
% 03/2025 also edited to store selected electrodes by subject for collab with George

clear; close all; clc

db = CCDTdatabase;
% db = CCDTdatabase_old;
Nsubj = height(db);

saveFile = 1;
savedir = '/mnt/sdb1/CCDT/CCDTscripts/CCDTnodeAnalyze dat/';
% savefname = 'CCDTnodeAnalyze_PLV_MH.mat'; % 'CCDTnodeAnalyze_PLI_MH.mat' 'CCDTnodeAnalyze_WPLI_MH.mat'

doZ = 0;
zType = [3]; % 1 = trial, 2 = electrode, 3 = freq, 'all' = whole matrix

connType = 0; % 0 = PLV, 1 = PLI, 2 = WPLI
if connType == 0
    savefname = 'CCDTnodeAnalyze_PLV_';
elseif connType == 1
    savefname = 'CCDTnodeAnalyze_PLI_';
elseif connType == 2
    savefname = 'CCDTnodeAnalyze_WPLI_';
end

use100 = 1
if use100
%     savefname = [savefname '100_pCAR_qCAR_hg70110_090924.mat'];    
    savefname = [savefname '100_pCAR_qCAR_hg70110_intra_090924.mat'];  
    savetablenameP = ['pow_selected_feats_preparatory_100dat.csv'];
    savetablenameQ = ['comm_selected_feats_preparatory_100dat.csv'];

else
    savefname = [savefname '80_pCAR_plvCAR_hg70110_30_031125.mat'];
    savetablenameP = ['pow_selected_feats_anticipatory.csv'];
    savetablenameQ = ['comm_selected_feats_anticipatory.csv'];

    percThresh = .3; % threshold w/ highest mean AUC from SVM analysis
end

%just nCom and pow with PLV, but has common median ref
% ddir = '/home/sharedMMRdata/CCDT Data/CCDT/CCDTanalysis_wCAR_noAmpOut_BR_072823/';
%pow and nCom for PLV, PLI, and wPLI, only has BP and CAR ref
% ddir = '/home/sharedMMRdata/CCDT Data/CCDT/CCDTanalysis_PLI_BR_081823/'; %brandon's data directory
% ddir = '/home/sharedMMRdata/CCDT Data/CCDT/CCDTanalysis_PLI_MH_082823/';
% rtdir = '/mnt/sdb1/CCDT/CCDT Data/CCDT/procData/'; % vivek's data directory
rtdir = '/mnt/sdb1/CCDT/orig_procData/'; % vivek's data directory
rtfnm = 'graphRT_sIall_500mspre'; %use graph metric versus PLV metric defined here
load([rtdir rtfnm '.mat'],'behavStruct'); %get RT and fast/slow trial info from this

pdir = '/mnt/sdb1/CCDT/procData072023/';
% pfname = 'powRT_101123_Bipolar_zAll'; %use for main pretrail analysis - 100% test fall data BP
% gfname = 'graphRT_090923_CAR_intra_MH';
% pfname = 'powRT_071023_Bipolar_BR_zAll';
% pfname = 'powRT_082823_CAR_BR_zAll'; % 100% test fall data CAR
% pfname = 'powRT_090224_BP_hg70150_1'; % reverted changes, new scripts
pfname = 'powRT_090924_CAR_hg70110.mat'; % old scripts
% pfname = 'powRT_090924_CAR_intra_hg70110.mat'; % old scripts intra
% pfname = 'powRT_090923_BP_intra_MH_zAll'; %intratrial
% pfname = 'powRT_031324_CAR_hg70-110'; %diff hg range pre
% pfname = 'powRT_050624_CAR_intra_hg70-110'; %diff hg range intra
% pfname = 'powRT_103123_Bipolar_zAll_null2';
% pfname = 'powRT_090923_CAR_intra_MH_zAll';
% gfname = 'graphRTpli_083123_CAR_BR';
% gfname = 'graphRT_071823_CAR_BR'; %use for main pretrial analysis - 100% test fall data
% gfname = 'graphRT_090224_CAR_pre_hg70150_1'; %reverted changes, new scripts
% gfname = 'graphRT_031025_CAR_PLV_hg70110'; %old scripts plv
gfname = 'graphRT_090324_CAR_hg70110'; %old scripts
% gfname = 'graphRT_090324_CAR_intra_hg70110'; %old scripts intra
% gfname = 'graphRT_102423_CAR_intra'; %intratrial
% gfname = 'graphRT_031024_CAR_hg70-150'; %diff hg range pre
% gfname = 'graphRT_050624_CAR_intra_hg70-110'; %diff hg range intra
% gfname = 'graphRT_103123_CAR_null2'; 

load([pdir pfname],"NFstruct","LTpow"); 
% load([rtdir 'powRT_sIall_500mspre'],"NFstruct","LTpow"); 
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","LTqexp");
% load([rtdir 'graphRT_sIall_500mspre'],"NFstruct","LTqexp");
NFstructQ = NFstruct;
clear NFstruct;

% load selected/significant features (only includes PLV data so far)
sigdir = '/mnt/sdb1/CCDT/CCDTscripts/80-20 dat/';
% fname8020 = '80-20-dat-pBP-qCAR-null2-103123.mat'; 
% fname8020 = '80-20-dat-pBP-qCAR-noZ-nullRT-101923.mat';  % pretrial
% fname8020 = '80-20-dat-pBP-qCAR-intra-102423.mat'; %intratrial
% fname8020 = '80-20-dat-pCAR-plvCAR-hg70110-031025.mat'; %pre diff hg range plv
% fname8020 = '80-20-dat-pCAR-qCAR-hg70110-090924.mat'; %pre diff hg range
fname8020 = '80-20-dat-pCAR-qCAR-hg70110-intra-090924.mat'; %intra diff hg range
if use100
    % get info from LTpow and LTqexp
    for isubj = 1:Nsubj
        sigNodesPow{isubj,1} = LTpow{isubj}(:,:,4)<=.05;
        sigNodesCom{isubj,1} = LTqexp{isubj}(:,:,4)<=.05;
    end
else
    load([sigdir fname8020],'sigPow','sigCom')
%     load([sigdir fname8020],'sigPowNdat','sigComNdat')
%     sigCom = sigComNdat;
%     sigPow = sigPowNdat;
    for isubj = 1:Nsubj
        sigNodesPow{isubj,1} = (sum(sigPow{isubj,1},3)/1000 > percThresh);
        sigNodesCom{isubj,1} = (sum(sigCom{isubj,1},3)/1000 > percThresh);
    end
end

%% analysis section
fPv = []; sPv = [];
fIDp = []; ssIDp = [];
fQv = []; sQv = [];
aPv = []; aQv = [];
fIDq = []; ssIDq = [];
gchP = []; gchQ = [];
for isubj = 1:Nsubj
    disp(['Subject: ' num2str(isubj)])

    % get iF, iS
    RT = behavStruct(isubj).vRT;
    iF = behavStruct(isubj).ifast;
    iS = behavStruct(isubj).islow;

    for j = 1:4 %fbands
        nComDat(:,:,j) = NFstructQ(isubj).fbands(j).NFqexp_nodal;
        nPowDat(:,:,j) = NFstructP(isubj).fbands(j).NFpow;
    end

    sigFeatP = sigNodesPow{isubj}; %nElec x nFreq
    nElecP = size(nPowDat,1);
    nTrl = size(nPowDat,2);
    nFreq = size(nPowDat,3);
    sigFeatQ = sigNodesCom{isubj};
    nElecCAR = size(nComDat,1);
    if doZ 
        nComDat = zscore(nComDat,[],zType);
        nPowDat = zscore(nPowDat,[],zType);
    end
    

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
%     a = 1;
    clear nComDat nPowDat
end
ssIDandDateQ= db(ssIDq,:);
ssIDandDateP = db(ssIDp,:);


if saveFile 
%     save([savedir savefname],'fPv','sPv','fIDp','ssIDp','gchP','fQv','sQv','aQv','aPv','fIDq','ssIDq','gchQ')
    selectedFeatureTable = table(ssIDp,ssIDandDateP(:,1),ssIDandDateP(:,2),gchP,fIDp,...
        'VariableNames',{'SessionNum','SubjectID','SessionDate','ChannelName','FrequencyBand'});
    writetable(selectedFeatureTable,[savedir savetablenameP])
    clear selectedFeatureTable
    selectedFeatureTable = table(ssIDq,ssIDandDateQ(:,1),ssIDandDateQ(:,2),gchQ,fIDq,...
        'VariableNames',{'SessionNum','SubjectID','SessionDate','ChannelName','FrequencyBand'});
    writetable(selectedFeatureTable,[savedir savetablenameQ])
end