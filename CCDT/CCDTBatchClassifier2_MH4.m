% function CCDTBatchClassifier2
% Train SVM classifier on combined neural data -- WM vs. GM with full
% network
%   VB 03/2020
%   MH 08/2023 - edited by MH to use data processed with original scripts
%   updated to work for null comparisons done with 80/20 data
%   use MH4 to loop through all individual regions.
clear; close all; clc
BTRChange = 27; %was 27

warning off
% parameters
pdir = '/mnt/sdb1/CCDT/procData072023/';
% pfname = 'powRT_090923_CAR_intra_MH_zAll';
% pfname = 'powRT_071023_Bipolar_BR_zAll';
% pfname = 'powRT_101123_Bipolar_zAll'; %main data to use
% pfname = 'powRT_031324_CAR_hg70-110';
% pfname = 'powRT_050624_CAR_intra_hg70-110'; %diff hg range CAR
pfname = 'powRT_090924_CAR_hg70110'; %from old scripts
% pfname = 'powRT_050624_BP_intra_hg70-110_removeOutl'; %diff hg range BP
% pfname = 'powRT_090923_BP_intra_MH_zAll';
% pfname = 'powRT_103123_Bipolar_zAll_null2';
% pfname = 'powRT_082823_CAR_BR_zAll';
% gfname = 'graphRT_071823_CAR_BR'; %main data to use
% gfname = 'graphRT_102423_CAR_intra';
% gfname = 'graphRT_103123_CAR_null2';
% gfname = 'graphRT_020424_CAR_PLV'; %plv instead of qexp
% gfname = 'graphRT_031024_CAR_hg70-110';
% gfname = 'graphRT_050624_CAR_intra_hg70-110'; %diff hg range
gfname = 'graphRT_090324_CAR_hg70110'; %from old scripts
% gfname = 'graphRTpli_083123_CAR_BR';
% gfname = 'graphRT_090923_CAR_intra_MH';
rtdir = '/mnt/sdb1/CCDT/CCDT Data/CCDT/procData/';
rtfname = 'graphRT_sIall_500mspre';
% pfnm = 'graphRT_sIall_500mspre'; %use graph metric versus PLV metric defined here
% powfnm = 'powRT_sIall_500mspre';
dir8020 = '/mnt/sdb1/CCDT/CCDTscripts/80-20 dat/';
% fname8020 = '80-20-dat-pBP-qCAR-noZ-nullRT-101923.mat'; % use for pre-trial and nullRT
% fname8020 = '80-20-dat-pBP-qCAR-intra-102423'; %intratrial
% fname8020 = '80-20-dat-pBP-qCAR-PLV-020424'; %plv
% fname8020 = '80-20-dat-pCAR-qCAR-hg70110-031324'; %pretrial diff hg range
% fname8020 = '80-20-dat-pBP-qCAR-hg70110-intra-050624'; %intratrial diff hg range
% fname8020 = '80-20-dat-pCAR-qCAR-hg70110-null-061324'; %shuffled RT null
fname8020 = '80-20-dat-pCAR-qCAR-hg70110-090924'; %pretrial from old scripts

plotSvNm = 'AUC80_thresh_CAR_hg70110_090924.mat';
saveForPlot = 1; %threshold selection plot/null comparison plots
Qpolarity = 1; %1 = CAR, 2 = BP
Ppolarity = 1; %1 = CAR, 2 = BP
useHalf = 0; % 0 = thirds vs. 1 = half
fbands = [3 12; 12 30; 30 55; 70 110]; %[3 12; 12 30; 30 55; 70 150];
SVMdat = 4; % 0 = 8fs (normal); 1 = 4fsqexp (or 4fsPLV); 2 = 4fspow; 3 = anatomically guided; 4 = anatomical nodes
% nodeType = 7; %if SVMdat = 4 --> 1 FrGM; 2 TmpGM; 3 ParGM; 4 TCWM; 5 FAWM; 6 TAWM; 7 ParWM; 8 ComWM
saveon = 0;
kfold=5;
minNodeType = 1; %0 for min/max slope, 1 for single min p value, 2 for multivariate regression
nullType = 0; % 0 for shuffle labels, 1 for simulated vRT, 2 for permuted vRT, 3 for random alteration of data
usePrevTrial = 0; % if nullType == 0, use previous trial (1) or randomly shuffle (0) for post extraction null
cvpart = 0;%[0.05:0.05:0.85]; %0 for use full set with kfold, otherwise specify fraction for temporal partition
setNewPthresh = 1; %1 = set a new P thresh instead of the file loaded
sigPThresh = 0.05; % if 1 above, this would be the new P threshold for feature selection
mvThresh = 0.1; % if using mv regression, set threshold for significant nodes
db = CCDTdatabase;
use100 = 0;
% thresh = .40;
if ~use100
    thresh = (0:5:100)/100; %set range of thresholds
    thresh_null = thresh;
%     thresh = .4; %set one threshold after selecting from range
%     thresh_null = .4; %nullType = 3: .15; nullType = 2: .4; nullType = 0: same as thresh
end

% load processed data

%loads: 
% behavStruct (contains RTs and iF, iS info) - still load this!
% LTmod (modularity data for all subjs)
% LTqexp (com data for all subjs)
% LTstr (strength data for all subjs)
% NFstruct (has a loooot of things in it - not sure how much needed)
    % fbands
    % gch - indices of good channels
    % gchlbl - labels of good channels
% sP (??)
% sPval (??)
% SVMstruct (results of SVM analysis?)
%     load([pdir pfnm '.mat']);
%     Nsubj = length(NFstruct);

for r = 1:8 %regions
    nodeType = r;
load([pdir pfname],"NFstruct","LTpow"); %
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","LTqexp");
NFstructQ = NFstruct;
clear NFstruct;
%load behavStruct RT data
load([rtdir rtfname],"behavStruct");
Nsubj = height(db);
load('all_WM.mat')
load patient_loc_120623

for iThresh = 1:length(thresh)
    percThresh = thresh(iThresh);
    percThresh_null = thresh_null(iThresh);
    disp(['Threshold: ' num2str(thresh(iThresh))])
    disp(['Null Threshold: ' num2str(thresh_null(iThresh))])
% load significant nodes from brandon's data
if use100
%     load([dir8020 '100-dat.mat']);
    % get info from LTpow and LTqexp
    for isubj = 1:Nsubj
        sigNodesPow{isubj,1} = LTpow{isubj}(:,:,4)<=.05;
        sigNodesCom{isubj,1} = LTqexp{isubj}(:,:,4)<=.05;
    end
else
    % for normal data load
    load([dir8020 fname8020],'sigPow','sigCom')  

    % for shuffled data null
%     load([dir8020 fname8020],'sigPowNdat','sigComNdat')
%     sigCom = sigComNdat;
%     sigPow = sigPowNdat;

    % for shuffled RT null
%     load([dir8020 fname8020],'sigPowNRT','sigComNRT')
%     sigCom = sigComNRT;
%     sigPow = sigPowNRT;

    for isubj = 1:Nsubj
        sigNodesPow{isubj,1} = (sum(sigPow{isubj,1},3)/1000 > percThresh);
        sigNodesCom{isubj,1} = (sum(sigCom{isubj,1},3)/1000 > percThresh);
    end
    if nullType == 2
        load([dir8020 fname8020],'nullShuffleRT','sigPowNRT','sigComNRT')
        for isubj = 1:Nsubj
            sigNodesPowNRT{isubj,1} = (sum(sigPowNRT{isubj,1},3)/1000 > percThresh_null);
            sigNodesComNRT{isubj,1} = (sum(sigComNRT{isubj,1},3)/1000 > percThresh_null);
        end 
    end
    if nullType == 3
        %idk when next time this option will get used again, but i
        %commented this stuff out so it should throw an error here so i 
        %remember to load in the right variables

%         load([dir8020 fname8020],'nullRandData','sigPowNdat','sigComNdat')
%         load([dir8020 fname8020],'nullRandData','sigPowNdatAdd','sigComNdatAdd')
        for isubj = 1:Nsubj
%             sigNodesPowNdat{isubj,1} = (sum(sigPowNdat{isubj,1},3)/1000 > percThresh_null);
%             sigNodesComNdat{isubj,1} = (sum(sigComNdat{isubj,1},3)/1000 > percThresh_null);
            sigNodesPowNdat{isubj,1} = (sum(sigPowNdatAdd{isubj,1},3)/1000 > percThresh_null);
            sigNodesComNdat{isubj,1} = (sum(sigComNdatAdd{isubj,1},3)/1000 > percThresh_null);
        end 
    end
    clear sigPow* sigCom*
end


LTnullStructRT20 = struct;
    
%     disp(['Iteration: ' num2str(jx)])
    for isubj = 1:Nsubj        
        disp(['Subject: ' num2str(isubj)])
        qexpNodes = cell(size(fbands,1),3); %raw channel, channel index, channel label
        powNodes = cell(size(fbands,1),3);
        cbehavStruct = behavStruct(isubj);
        cNFstructP = NFstructP(isubj);
        cNFstructQ = NFstructQ(isubj);
        cCHpow = cNFstructP.gch;
        cCH = cNFstructQ.gch;
        cgchlblpow = cNFstructP.gchlbl;
        cgchlbl = cNFstructQ.gchlbl;

        cRT = cbehavStruct.vRT; % reaction time (ms)
        goodTrl = ~(cRT<50 | cRT>999);

        
        %Behavior setup
        if ~useHalf
            iF = cbehavStruct.ifast;
            iS = cbehavStruct.islow;
%             iF = cbehavStruct.ifast(goodTrl);
%             iS = cbehavStruct.islow(goodTrl);
        else
            cTc = (cRT>50&cRT<999);
            [~,b] = sort(cRT, 'ascend');
            cB = (b(cTc(b)));
            iF = zeros(length(cTc),1);
            iS = zeros(length(cTc),1);
            iF(cB(1:round(sum(cTc)/2))) = 1;
            iS(cB(end-round(sum(cTc)/2)+1:end)) = 1;
            iF=logical(iF);
            iS=logical(iS);
        end
        cN = cell(length(cRT),1);
        cN(iF) = {'Fast'};
        cN(iS) = {'Slow'};
        cN(~iF&~iS) = [];
        cNRT = cRT(iF|iS);
               
        if nullType==1
            nullRT = rand(length(cRT),1).*999;
            inr = (cRT<50 | cRT>999);
        elseif nullType==2
%                 nullRT = cRT(randperm(length(cRT)));
            nullRT = nullShuffleRT{isubj};
            inr = (nullRT<50 | nullRT>999);
        elseif nullType==3
            nullRT = cRT;
            inr = (nullRT<50 | nullRT>999); %bad RT
        else %nullType == 0 shuffle label null
            nullRT = [];
            if ~usePrevTrial
                ncN = cN(randperm(length(cN)));
            else
                nullRT(1) = cRT(end);
                for i=2:length(cRT)
                    nullRT(i) = cRT(i-1);
                end
            end
        end

        if nullType > 0 || (~nullType && usePrevTrial)
            if nullType>0
                ncTc = (nullRT>50&nullRT<999);
            else
                ncTc = (nullRT>0&nullRT<999); 
            end
            [~,b] = sort(nullRT, 'ascend');
            cB = (b(ncTc(b)));
            iFn = zeros(length(ncTc),1);
            iSn = zeros(length(ncTc),1);
            if ~useHalf
                iFn(cB(1:round(sum(ncTc)/3))) = 1;
                iSn(cB(end-round(sum(ncTc)/3)+1:end)) = 1;
            else
                iFn(cB(1:round(sum(ncTc)/2))) = 1;
                iSn(cB(end-round(sum(ncTc)/2)+1:end)) = 1;
            end
            iFn=logical(iFn);
            iSn=logical(iSn);
            ncN = cell(length(nullRT),1);
            ncN(iFn) = {'Fast'};
            ncN(iSn) = {'Slow'};
            ncN(~iFn&~iSn) = [];
        end
    
        
        tNFpow = cNFstructP.fbands(1).NFpow;
        bNFpow = cNFstructP.fbands(2).NFpow;
        lgNFpow = cNFstructP.fbands(3).NFpow;
        hgNFpow = cNFstructP.fbands(4).NFpow;
        
        tNFqexp = cNFstructQ.fbands(1).NFqexp_nodal;
        bNFqexp = cNFstructQ.fbands(2).NFqexp_nodal;
        lgNFqexp = cNFstructQ.fbands(3).NFqexp_nodal;
        hgNFqexp = cNFstructQ.fbands(4).NFqexp_nodal;
         
%         sPpow = powDat.sP;
%         LTpow = powDat.LTpow;
        [~,tPmSi] = min(LTpow{isubj}(:,1,1)); [~,tPmaSi] = max(LTpow{isubj}(:,1,1));
        [~,bPmSi] = min(LTpow{isubj}(:,2,1)); [~,bPmaSi] = max(LTpow{isubj}(:,2,1));
        [~,lgPmSi] = min(LTpow{isubj}(:,3,1)); [~,lgPmaSi] = max(LTpow{isubj}(:,3,1));
        [~,hgPmSi] = min(LTpow{isubj}(:,4,1)); [~,hgPmaSi] = max(LTpow{isubj}(:,4,1));
        
        [~,tQmSi] = min(LTqexp{isubj}(:,1,1)); [~,tQmaSi] = max(LTqexp{isubj}(:,1,1));
        [~,bQmSi] = min(LTqexp{isubj}(:,2,1)); [~,bQmaSi] = max(LTqexp{isubj}(:,2,1));
        [~,lgQmSi] = min(LTqexp{isubj}(:,3,1)); [~,lgQmaSi] = max(LTqexp{isubj}(:,3,1));
        [~,hgQmSi] = min(LTqexp{isubj}(:,4,1)); [~,hgQmaSi] = max(LTqexp{isubj}(:,4,1));

        
        
        if setNewPthresh
            clear useQnodes usePnodes 
            for i=1:length(fbands)
                %use 100 feature selection
%                     useQnodes{i} = find(LTqexp{isubj}(:,i,4)<sigPThresh);
%                     usePnodes{i} = find(LTpow{isubj}(:,i,4)<sigPThresh);
%                     nsQnodes{i} = find(LTqexp{isubj}(:,i,4)>sigPThresh);
%                     nsPnodes{i} = find(LTpow{isubj}(:,i,4)>sigPThresh);

                %use 80/20 feature selction
                if SVMdat == 4 
                    %just get TCWM nodes
                    nodes_in_reg = patient_loc(Qpolarity).session(isubj).gm_wm_rois == nodeType;
                    useQnodes{i} = find(sigNodesCom{isubj}(:,i) & nodes_in_reg);
                    nsQnodes{i} = find(~sigNodesCom{isubj}(:,i) | ~nodes_in_reg);

                    nodes_in_reg = patient_loc(Ppolarity).session(isubj).gm_wm_rois == nodeType;                    
                    usePnodes{i} = find(sigNodesPow{isubj}(:,i) & nodes_in_reg);
                    nsPnodes{i} = find(~sigNodesPow{isubj}(:,i) | ~nodes_in_reg);
                else
                    useQnodes{i} = find(sigNodesCom{isubj}(:,i));
                    usePnodes{i} = find(sigNodesPow{isubj}(:,i));
                    nsQnodes{i} = find(~sigNodesCom{isubj}(:,i));
                    nsPnodes{i} = find(~sigNodesPow{isubj}(:,i));
                end 
                if nullType==2
                    nullNodesQ{i} = find(sigNodesComNRT{isubj}(:,i));
                    nullNodesP{i} = find(sigNodesPowNRT{isubj}(:,i));
                end
                if nullType==3
                    nullNodesQ{i} = find(sigNodesComNdat{isubj}(:,i));
                    nullNodesP{i} = find(sigNodesPowNdat{isubj}(:,i));
                end
            end
        end
              
        %mv regression set up
        try
            [~, sMVt] = robustfit([tNFqexp(useQnodes{1},(iF|iS))'],cRT(iF|iS));
            mvQnodes{1} = useQnodes{1}(sMVt.p(2:end)<mvThresh);
        catch
            mvQnodes{1} = NaN;
        end
        try
            [~, sMVb] = robustfit([bNFqexp(useQnodes{2},(iF|iS))'],cRT(iF|iS));
            mvQnodes{2} = useQnodes{2}(sMVb.p(2:end)<mvThresh);
        catch
            mvQnodes{2} = NaN;
        end
        
        try
            [~, sMVlg] = robustfit([lgNFqexp(useQnodes{3},(iF|iS))'],cRT(iF|iS));
            mvQnodes{3} = useQnodes{3}(sMVlg.p(2:end)<mvThresh);
        catch
            mvQnodes{3} = NaN;
        end
        
        try
            [~, sMVhg] = robustfit([hgNFqexp(useQnodes{4},(iF|iS))'],cRT(iF|iS));
            mvQnodes{4} = useQnodes{4}(sMVhg.p(2:end)<mvThresh);
        catch
            mvQnodes{4} = NaN;
        end
        
        try
            [~, sMVtp] = robustfit([tNFpow(usePnodes{1},(iF|iS))'],cRT(iF|iS));
            mvPnodes{1} = usePnodes{1}(sMVtp.p(2:end)<mvThresh);
        catch
            mvPnodes{1} = NaN;
        end
        try
            [~, sMVbp] = robustfit([bNFpow(usePnodes{2},(iF|iS))'],cRT(iF|iS));
            mvPnodes{2} = usePnodes{2}(sMVbp.p(2:end)<mvThresh);
        catch
            mvPnodes{2} = NaN;
        end
        
        try
            [~, sMVlgp] = robustfit([lgNFpow(usePnodes{3},(iF|iS))'],cRT(iF|iS));
            mvPnodes{3} = usePnodes{3}(sMVlgp.p(2:end)<mvThresh);
        catch
            mvPnodes{3} = NaN;
        end
        try
            [~, sMVhgp] = robustfit([hgNFpow(usePnodes{4},(iF|iS))'],cRT(iF|iS));
            mvPnodes{4} = usePnodes{4}(sMVhgp.p(2:end)<mvThresh);
        catch
            mvPnodes{4} = NaN;
        end
        
        
        
        
        totalSigNodesQ_thresh(isubj,:) = cellfun(@length,useQnodes);
        totalSigNodesP_thresh(isubj,:) = cellfun(@length,usePnodes);
        totalSigNodesQmv_thresh(isubj,:) = cellfun(@length,mvQnodes);
        totalSigNodesPmv_thresh(isubj,:) = cellfun(@length,mvPnodes);
        if nullType>0
            totalSigNodesQnull_thresh(isubj,:) = cellfun(@length,nullNodesQ);
            totalSigNodesPnull_thresh(isubj,:) = cellfun(@length,nullNodesP);
        end
        
        if SVMdat==1
            dat = [tNFqexp(useQnodes{1},(iF|iS))' bNFqexp(useQnodes{2},(iF|iS))' lgNFqexp(useQnodes{3},(iF|iS))' hgNFqexp(useQnodes{4},(iF|iS))'];
            nsDat = [tNFqexp(nsQnodes{1},(iF|iS))' bNFqexp(nsQnodes{2},(iF|iS))' lgNFqexp(nsQnodes{3},(iF|iS))' hgNFqexp(nsQnodes{4},(iF|iS))'];
%             if minNodeType==1
%                 datMinNode = [tNFqexp((sP{isubj,1}(1)),(iF|iS))' bNFqexp((sP{isubj,2}(1)),(iF|iS))' lgNFqexp((sP{isubj,3}(1)),(iF|iS))' hgNFqexp((sP{isubj,4}(1)),(iF|iS))'];
%             elseif minNodeType==2
%                 datMinNode = [tNFqexp(mvQnodes{1},(iF|iS))' bNFqexp(mvQnodes{2},(iF|iS))' lgNFqexp(mvQnodes{3},(iF|iS))' hgNFqexp(mvQnodes{4},(iF|iS))'];
%             else
%                 datMinNode = [tNFqexp([tQmSi tQmaSi],(iF|iS))' bNFqexp([bQmSi bQmaSi],(iF|iS))' lgNFqexp([lgQmSi lgQmaSi],(iF|iS))' hgNFqexp([hgQmSi hgQmaSi],(iF|iS))'];
%             end
            
            
        elseif SVMdat==2
            dat = [tNFpow(usePnodes{1},(iF|iS))' bNFpow(usePnodes{2},(iF|iS))' lgNFpow(usePnodes{3},(iF|iS))' hgNFpow(usePnodes{4},(iF|iS))'];
            nsDat = [tNFpow(nsPnodes{1},(iF|iS))' bNFpow(nsPnodes{2},(iF|iS))' lgNFpow(nsPnodes{3},(iF|iS))' hgNFpow(nsPnodes{4},(iF|iS))'];
%             if minNodeType==1
%                 datMinNode = [tNFpow((sPpow{isubj}(1,1)),(iF|iS))' bNFpow((sPpow{isubj}(2,1)),(iF|iS))' lgNFpow((sPpow{isubj}(3,1)),(iF|iS))' hgNFpow((sPpow{isubj}(4,1)),(iF|iS))' ];
%             elseif minNodeType==2
%                 datMinNode = [tNFpow(mvPnodes{1},(iF|iS))' bNFpow(mvPnodes{2},(iF|iS))' lgNFpow(mvPnodes{3},(iF|iS))' hgNFpow(mvPnodes{4},(iF|iS))'];
%             else
%                 datMinNode = [tNFpow([tPmSi tPmaSi],(iF|iS))' bNFpow([bPmSi bPmaSi],(iF|iS))' lgNFpow([lgPmSi lgPmaSi],(iF|iS))' hgNFpow([hgPmSi hgPmaSi],(iF|iS))'];
%             end
            
            
        elseif SVMdat==3
            load('patient_loc_101421.mat')
            % Based on heat map
            usePnodes{1}=[];
            usePnodes{2} = (patient_loc(2).session(isubj).wm==3);
            usePnodes{3} = [];
            usePnodes{4} = [];
            useQnodes{1} = (patient_loc(1).session(isubj).wm==1);
            useQnodes{2} = (patient_loc(1).session(isubj).wm==2|patient_loc(1).session(isubj).wm==3|patient_loc(1).session(isubj).seg==3);
            useQnodes{3} = [];
            useQnodes{4} = (patient_loc(1).session(isubj).seg==1);
            dat = [tNFpow(usePnodes{1},(iF|iS))' bNFpow(usePnodes{2},(iF|iS))' lgNFpow(usePnodes{3},(iF|iS))' hgNFpow(usePnodes{4},(iF|iS))' tNFqexp(useQnodes{1},(iF|iS))' bNFqexp(useQnodes{2},(iF|iS))' lgNFqexp(useQnodes{3},(iF|iS))' hgNFqexp(useQnodes{4},(iF|iS))'];
            nsDat = [tNFpow(nsPnodes{1},(iF|iS))' bNFpow(nsPnodes{2},(iF|iS))' lgNFpow(nsPnodes{3},(iF|iS))' hgNFpow(nsPnodes{4},(iF|iS))' tNFqexp(nsQnodes{1},(iF|iS))' bNFqexp(nsQnodes{2},(iF|iS))' lgNFqexp(nsQnodes{3},(iF|iS))' hgNFqexp(nsQnodes{4},(iF|iS))'];
%             if minNodeType==1
%                 datMinNode = [tNFpow((sPpow{isubj}(1,1)),(iF|iS))' bNFpow((sPpow{isubj}(2,1)),(iF|iS))' lgNFpow((sPpow{isubj}(3,1)),(iF|iS))' hgNFpow((sPpow{isubj}(4,1)),(iF|iS))' tNFqexp((sP{isubj,1}(1)),(iF|iS))' bNFqexp((sP{isubj,2}(1)),(iF|iS))' lgNFqexp((sP{isubj,3}(1)),(iF|iS))' hgNFqexp((sP{isubj,4}(1)),(iF|iS))'];
%             elseif minNodeType==2
%                 datMinNode = [tNFpow(mvPnodes{1},(iF|iS))' bNFpow(mvPnodes{2},(iF|iS))' lgNFpow(mvPnodes{3},(iF|iS))' hgNFpow(mvPnodes{4},(iF|iS))' tNFqexp(mvQnodes{1},(iF|iS))' bNFqexp(mvQnodes{2},(iF|iS))' lgNFqexp(mvQnodes{3},(iF|iS))' hgNFqexp(mvQnodes{4},(iF|iS))'];
%             else
%                 datMinNode = [tNFpow([tPmSi tPmaSi],(iF|iS))' bNFpow([bPmSi bPmaSi],(iF|iS))' lgNFpow([lgPmSi lgPmaSi],(iF|iS))' hgNFpow([hgPmSi hgPmaSi],(iF|iS))' tNFqexp([tQmSi tQmaSi],(iF|iS))' bNFqexp([bQmSi bQmaSi],(iF|iS))' lgNFqexp([lgQmSi lgQmaSi],(iF|iS))' hgNFqexp([hgQmSi hgQmaSi],(iF|iS))'];
%             end
            
        else
            dat = [tNFpow(usePnodes{1},(iF|iS))' bNFpow(usePnodes{2},(iF|iS))' lgNFpow(usePnodes{3},(iF|iS))' hgNFpow(usePnodes{4},(iF|iS))' tNFqexp(useQnodes{1},(iF|iS))' bNFqexp(useQnodes{2},(iF|iS))' lgNFqexp(useQnodes{3},(iF|iS))' hgNFqexp(useQnodes{4},(iF|iS))'];
            nsDat = [tNFpow(nsPnodes{1},(iF|iS))' bNFpow(nsPnodes{2},(iF|iS))' lgNFpow(nsPnodes{3},(iF|iS))' hgNFpow(nsPnodes{4},(iF|iS))' tNFqexp(nsQnodes{1},(iF|iS))' bNFqexp(nsQnodes{2},(iF|iS))' lgNFqexp(nsQnodes{3},(iF|iS))' hgNFqexp(nsQnodes{4},(iF|iS))'];
%             if minNodeType==1
%                 datMinNode = [tNFpow((sPpow{isubj}(1,1)),(iF|iS))' bNFpow((sPpow{isubj}(2,1)),(iF|iS))' lgNFpow((sPpow{isubj}(3,1)),(iF|iS))' hgNFpow((sPpow{isubj}(4,1)),(iF|iS))' tNFqexp((sP{isubj,1}(1)),(iF|iS))' bNFqexp((sP{isubj,2}(1)),(iF|iS))' lgNFqexp((sP{isubj,3}(1)),(iF|iS))' hgNFqexp((sP{isubj,4}(1)),(iF|iS))'];
%             elseif minNodeType==2
%                 datMinNode = [tNFpow(mvPnodes{1},(iF|iS))' bNFpow(mvPnodes{2},(iF|iS))' lgNFpow(mvPnodes{3},(iF|iS))' hgNFpow(mvPnodes{4},(iF|iS))' tNFqexp(mvQnodes{1},(iF|iS))' bNFqexp(mvQnodes{2},(iF|iS))' lgNFqexp(mvQnodes{3},(iF|iS))' hgNFqexp(mvQnodes{4},(iF|iS))'];
%             else
%                 datMinNode = [tNFpow([tPmSi tPmaSi],(iF|iS))' bNFpow([bPmSi bPmaSi],(iF|iS))' lgNFpow([lgPmSi lgPmaSi],(iF|iS))' hgNFpow([hgPmSi hgPmaSi],(iF|iS))' tNFqexp([tQmSi tQmaSi],(iF|iS))' bNFqexp([bQmSi bQmaSi],(iF|iS))' lgNFqexp([lgQmSi lgQmaSi],(iF|iS))' hgNFqexp([hgQmSi hgQmaSi],(iF|iS))'];
%             end
            
            
        end
        
        datAllNodes = [tNFpow(:,(iF|iS))' bNFpow(:,(iF|iS))' lgNFpow(:,(iF|iS))' hgNFpow(:,(iF|iS))' tNFqexp(:,(iF|iS))' bNFqexp(:,(iF|iS))' lgNFqexp(:,(iF|iS))' hgNFqexp(:,(iF|iS))'];
        
        
        if nullType>0
            if SVMdat==1
                nullDat = [tNFqexp(nullNodesQ{1},(iFn|iSn))' bNFqexp(nullNodesQ{2},(iFn|iSn))' lgNFqexp(nullNodesQ{3},(iFn|iSn))' hgNFqexp(nullNodesQ{4},(iFn|iSn))'];
            elseif SVMdat==2
                nullDat = [tNFpow(nullNodesP{1},(iFn|iSn))' bNFpow(nullNodesP{2},(iFn|iSn))' lgNFpow(nullNodesP{3},(iFn|iSn))' hgNFpow(nullNodesP{4},(iFn|iSn))'];
            else
                if nullType==3
                    % this should be set up to use noisy to data.
                    % change out the commented sections to subtract data
                    % off to just look at noise or toggle multiplied/added noise

                    cnullRandDataPow = nullRandData{isubj,1}; %multiplied noise 
                    cnullRandDataCom = nullRandData{isubj,2}; %multiplied noise
%                     cnullRandDataPow = nullRandData{isubj,3}; %added noise
%                     cnullRandDataCom = nullRandData{isubj,4}; %added noise

                    tNFpow = cnullRandDataPow(:,:,1);
                    bNFpow = cnullRandDataPow(:,:,2);
                    lgNFpow = cnullRandDataPow(:,:,3);
                    hgNFpow = cnullRandDataPow(:,:,4);
                    
                    tNFqexp = cnullRandDataCom(:,:,1);
                    bNFqexp = cnullRandDataCom(:,:,2);
                    lgNFqexp = cnullRandDataCom(:,:,3);
                    hgNFqexp = cnullRandDataCom(:,:,4);

                    %subtract data to just get the noise added
%                     tNFpow = cnullRandDataPow(:,:,1) - cNFstructP.fbands(1).NFpow(:,goodTrl);
%                     bNFpow = cnullRandDataPow(:,:,2) - cNFstructP.fbands(2).NFpow(:,goodTrl);
%                     lgNFpow = cnullRandDataPow(:,:,3) - cNFstructP.fbands(3).NFpow(:,goodTrl);
%                     hgNFpow = cnullRandDataPow(:,:,4) - cNFstructP.fbands(4).NFpow(:,goodTrl);
%                     
%                     tNFqexp = cnullRandDataCom(:,:,1) - cNFstructQ.fbands(1).NFqexp_nodal(:,goodTrl);
%                     bNFqexp = cnullRandDataCom(:,:,2) - cNFstructQ.fbands(2).NFqexp_nodal(:,goodTrl);
%                     lgNFqexp = cnullRandDataCom(:,:,3) - cNFstructQ.fbands(3).NFqexp_nodal(:,goodTrl);
%                     hgNFqexp = cnullRandDataCom(:,:,4) - cNFstructQ.fbands(4).NFqexp_nodal(:,goodTrl);


                    iFn = iFn(goodTrl);
                    iSn = iSn(goodTrl);
                end
                
                nullDat = [tNFpow(nullNodesP{1},(iFn|iSn))' bNFpow(nullNodesP{2},(iFn|iSn))' lgNFpow(nullNodesP{3},(iFn|iSn))' hgNFpow(nullNodesP{4},(iFn|iSn))' tNFqexp(nullNodesQ{1},(iFn|iSn))' bNFqexp(nullNodesQ{2},(iFn|iSn))' lgNFqexp(nullNodesQ{3},(iFn|iSn))' hgNFqexp(nullNodesQ{4},(iFn|iSn))'];
            end
        end
        
        trainSet = zeros(length(cN),1);
        for iij = 1:length(cvpart)
            ccvpart = cvpart(iij);
            if ccvpart>0
                trainSet(1:length(cN)*ccvpart) = 1;
            else
                try
                    mdlSVMcv = fitcsvm(dat,cN,'KFold',kfold);
                    mdlSVMcvNS = fitcsvm(nsDat(:,randi(size(nsDat,2),size(dat,2))),cN,'KFold',kfold);
                    %     mdlSVMcv = fitclinear(dat,cN,'learner', 'svm','KFold',kfold,'regularization','ridge');
                catch
                    disp(['No significant feature nodes for subject session ' num2str(isubj) ' -- skipping...'])
                    break
                end
                [score_lbl,score_svm] = kfoldPredict(mdlSVMcv);
                %                 [X,Y,T,tAUC(ixi,:)] = perfcurve(cN,score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                [X,Y,T,AUC] = perfcurve(cN,score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                [score_lblNS,score_svmNS] = kfoldPredict(mdlSVMcvNS);
                %                 [Xns,Yns,Tns,tAUCns(ixi,:)] = perfcurve(cN,score_svmNS(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                [Xns,Yns,Tns,AUCns] = perfcurve(cN,score_svmNS(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                %             AUC = mean(tAUC);
                %             AUCns = mean(tAUCns);
                
                %logistic regression comparison
                mdlGLMcv = fitclinear(dat,cN,'learner','logistic','kfold',kfold,'regularization','ridge');
                [~,score_glm] = kfoldPredict(mdlGLMcv);
                [Xglm,Yglm,Tglm,AUCglm] = perfcurve(cN,score_glm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                mdlGLMcvNS = fitclinear(nsDat,cN,'learner','logistic','kfold',kfold,'regularization','ridge');
                [~,score_glmNS] = kfoldPredict(mdlGLMcvNS);
                [XglmNS,YglmNS,TglmNS,AUCglmNS] = perfcurve(cN,score_glmNS(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                
                
                % null comparison
                if ~nullType
                    %dat and ncN not same size now
                    mdlSVMcvnull = fitcsvm(dat,ncN,'KFold',kfold);
                    mdlGLMcvnull = fitclinear(dat,ncN,'learner','logistic','kfold',kfold, 'regularization','ridge');
                else
                    %nullDat and ncN are right size now
                    if ~isempty(nullDat)
                    mdlSVMcvnull = fitcsvm(nullDat,ncN,'KFold',kfold);
                    mdlGLMcvnull = fitclinear(nullDat,ncN,'learner','logistic','kfold',kfold, 'regularization','ridge');
                    end
                end
                if nullType >0
                    if ~isempty(nullDat)
                        [score_lblnull,score_svmnull] = kfoldPredict(mdlSVMcvnull);
                        [Xnull,Ynull,Tnull,AUCnull] = perfcurve(ncN,score_svmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                        
                        [~,score_glmnull] = kfoldPredict(mdlGLMcvnull);
                        [Xglmnull,Yglmnull,Tglmnull,AUCglmnull] = perfcurve(ncN,score_glmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                    else
                        score_lblnull = NaN; score_svmnull = NaN;
                        Xnull = NaN; Ynull = NaN; Tnull = NaN; AUCnull = NaN;
                        score_glmnull = NaN;
                        Xglmnull = NaN; Yglmnull = NaN; Tglmnull = NaN; AUCglmnull = NaN;
                    end
                else
                    [score_lblnull,score_svmnull] = kfoldPredict(mdlSVMcvnull);
                    [Xnull,Ynull,Tnull,AUCnull] = perfcurve(ncN,score_svmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                    
                    [~,score_glmnull] = kfoldPredict(mdlGLMcvnull);
                    [Xglmnull,Yglmnull,Tglmnull,AUCglmnull] = perfcurve(ncN,score_glmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                end
                %linear regression comparison
                mdlLM = fitrlinear(dat,cNRT,'learner','leastsquares','regularization','ridge','kfold',kfold);
                cNRTnull = cNRT(randperm(length(cNRT)));
                mdlLMnull = fitrlinear(dat,cNRTnull,'learner','leastsquares','regularization','ridge','kfold',kfold);
                %     mdlLMmin = fitlm(datMinNode,cNRT);
                %     mdlLMnullmin = fitlm(datMinNode,cNRTnull);
                %      mdlLM = fitlm(dat,cNRT);
                %     mdlLMnull = fitlm(dat,cNRTnull);
                
                
                % SVM with ridge and lasso
                mdlSVML2 = fitclinear(dat,cN,'learner','svm','regularization','ridge', 'kfold',kfold);
                [score_lblL2,score_svmL2] = kfoldPredict(mdlSVML2);
                [XL2,YL2,TL2,AUCL2] = perfcurve(cN,score_svmL2(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                mdlSVML1 = fitclinear(dat,cN,'learner','svm','regularization','lasso', 'kfold',kfold);
                [score_lblL2,score_svmL1] = kfoldPredict(mdlSVML1);
                [XL1,YL1,TL1,AUCL1] = perfcurve(cN,score_svmL1(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                %         mdlSVMoptimize = fitclinear(dat,cN,'OptimizeHyperparameters','all');
                
                % test all Nodes and single min node
                mdlSVMcv_allNodes = fitclinear(datAllNodes,cN,'learner', 'svm','KFold',kfold,'regularization','ridge');
                
                [~,score_svmAllNodes] = kfoldPredict(mdlSVMcv_allNodes);
                [Xall,Yall,Tall,AUCall] = perfcurve(cN,score_svmAllNodes(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
%                 try
%                     mdlSVMcv_minNode = fitclinear(datMinNode,cN,'learner', 'svm','KFold',kfold,'regularization','ridge');
%                 catch
%                     disp(['No significant multivariate nodes for subject session ' num2str(isubj) ' -- skipping...'])
%                     break
%                 end
%                 
%                 %     mdlSVMcv_minNode = fitcsvm(datMinNode,cN,'KFold',kfold);
%                 [~,score_svmMinNode] = kfoldPredict(mdlSVMcv_minNode);
%                 [Xmin,Ymin,Tmin,AUCmin] = perfcurve(cN,score_svmMinNode(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
            end
            
            trainSet = logical(trainSet);
            if ccvpart>0
                mdlSVMcv = fitcsvm(dat(trainSet,:),cN(trainSet));
                [score_lbl,score_svm] = predict(mdlSVMcv, dat(~trainSet,:));
                [X,Y,T,AUC(iij,:)] = perfcurve(cN(~trainSet),score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                if ~nullType
                    mdlSVMcvnull = fitcsvm(dat(trainSet,:),ncN(trainSet));
                else
                    mdlSVMcvnull = fitcsvm(nullDat(trainSet,:),ncN(trainSet));
                end
                [score_lblnull,score_svmnull] = predict(mdlSVMcvnull,dat(~trainSet,:));
                [Xnull,Ynull,Tnull,AUCnull(iij,:)] = perfcurve(ncN(~trainSet),score_svmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                mdlSVMcv_allNodes = fitcsvm(datAllNodes(trainSet,:),cN(trainSet));
                
                [~,score_svmAllNodes] = predict(mdlSVMcv_allNodes,datAllNodes(~trainSet,:));
                [Xall,Yall,Tall,AUCall(iij,:)] = perfcurve(cN(~trainSet),score_svmAllNodes(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                mdlSVMcv_minNode = fitcsvm(datMinNode(trainSet,:),cN(trainSet));
                [~,score_svmMinNode] = predict(mdlSVMcv_minNode, datMinNode(~trainSet,:));
                [Xmin,Ymin,Tmin,AUCmin(iij,:)] = perfcurve(cN(~trainSet),score_svmMinNode(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
            end
        end
        
        if exist("AUC")
        AUCasCVpart_inf(isubj,:,:) = AUC;
%         AUCasCVpart_min(isubj,:,:) = AUCmin;
        AUCasCVpart_all(isubj,:,:) = AUCall;
        AUCasCVpart_null(isubj,:,:) = AUCnull;
        AUCasCVpart_NS(isubj,:,:) = AUCns;
        
        
        batchSVMstruct(isubj).mdl.SVM = mdlSVMcv;
        
        %     batchSVMstruct(isubj).mdl.SVML1 = mdlSVML1;
        %     batchSVMstruct(isubj).mdl.SVML2 = mdlSVML2;
        %     batchSVMstruct(isubj).mdl.glm = mdlGLMcv;
        %     batchSVMstruct(isubj).mdl.lm = mdlLM;
        %     %         batchSVMstruct(isubj).mdl.optimize = mdlSVMoptimize;
        %     batchSVMstruct(isubj).mdl.SVMallNodes = mdlSVMcv_allNodes;
        %     batchSVMstruct(isubj).mdl.SVMminNode = mdlSVMcv_minNode;
        
        
        
        
        batchSVMstruct(isubj).qexpNodes = useQnodes;
        batchSVMstruct(isubj).powNodes = usePnodes;
        batchSVMstruct(isubj).dataIn = dat;
        batchSVMstruct(isubj).score_lbl = score_lbl;
        batchSVMstruct(isubj).score_svm = score_svm;
        batchSVMstruct(isubj).classNames = cN;
        batchSVMstruct(isubj).RT = cRT(iF|iS);
        batchSVMstruct(isubj).AUC.SVM = AUC;
        batchSVMstruct(isubj).X.SVM = X;
        batchSVMstruct(isubj).Y.SVM = Y;
        batchSVMstruct(isubj).T.SVM = T;
        
        
        
        batchSVMstruct(isubj).null.mdl = mdlSVMcvnull;
        batchSVMstruct(isubj).null.score_lbl = score_lblnull;
        batchSVMstruct(isubj).null.score_svm = score_svmnull;
        batchSVMstruct(isubj).null.classNames = ncN;
        batchSVMstruct(isubj).null.AUC = AUCnull;
        batchSVMstruct(isubj).null.X = Xnull;
        batchSVMstruct(isubj).null.Y = Ynull;
        batchSVMstruct(isubj).null.T = Tnull;
        
        if length(thresh)==1
            %figure 2D
            figure(1)
            set(gcf,'Color','white','Units','inches','Position',[2,2,5,3.5])
            hold on;
            % plot(X(:,1),Y(:,1), Y(:,2),Y(:,3),X(:,2),X(:,3),'b')
            % plot(Xnull(:,1),Ynull(:,1), Ynull(:,2),Ynull(:,3),Xnull(:,2),Xnull(:,3), 'k')
            % plot(X(:,1),Y(:,1),'b')
            % plot(Xnull(:,1),Ynull(:,1), 'k')
            shadedErrorBar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),'b');
            shadedErrorBar(Xnull(:,1),Ynull(:,1),Ynull(:,1)-Ynull(:,2),'k');
            plot([0,1],[0,1],'--','Color','black');
            hold off
            xlim([0,1]); xticks([0,.5,1]);
            ylim([0,1.1]); yticks([0,.5,1]);
            xlabel('True Positive Rate')
            ylabel('False Positive Rate')
            
            figure(2)
            set(gcf,'Color','white')
            sgtitle(['Behavioral (blue) vs Null (gray) ROC, nullType = ' num2str(nullType)])
            subplot(3,9,isubj)
            hold on
            shadedErrorBar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),'b');
            shadedErrorBar(Xnull(:,1),Ynull(:,1),Ynull(:,1)-Ynull(:,2),'k');
            title(['Subject ' num2str(isubj)])
            plot([0,1],[0,1],'r')
        end
        
        batchSVMstruct(isubj).AUC.SVMallNodes = AUCall;
        batchSVMstruct(isubj).X.SVMallNodes = Xall;
        batchSVMstruct(isubj).Y.SVMallNodes = Yall;
        batchSVMstruct(isubj).T.SVMallNodes = Tall;
        
%         batchSVMstruct(isubj).AUC.SVMminNode = AUCmin;
%         batchSVMstruct(isubj).X.SVMminNode = Xmin;
%         batchSVMstruct(isubj).Y.SVMminNode = Ymin;
%         batchSVMstruct(isubj).T.SVMminNode = Tmin;
        
        batchSVMstruct(isubj).AUC.SVMns = AUCns;
        batchSVMstruct(isubj).X.SVMns = Xns;
        batchSVMstruct(isubj).Y.SVMns = Yns;
        batchSVMstruct(isubj).T.SVMns = Tns;
        
        if ~cvpart
            batchSVMstruct(isubj).AUC.glm = AUCglm;
            batchSVMstruct(isubj).X.glm = Xglm;
            batchSVMstruct(isubj).Y.glm = Yglm;
            batchSVMstruct(isubj).T.glm = Tglm;
            
            batchSVMstruct(isubj).AUC.SVML2 = AUCL2;
            batchSVMstruct(isubj).X.SVML2 = XL2;
            batchSVMstruct(isubj).Y.SVML2 = YL2;
            batchSVMstruct(isubj).T.SVML2 = TL2;
            
            batchSVMstruct(isubj).AUC.SVML1 = AUCL1;
            batchSVMstruct(isubj).X.SVML1 = XL1;
            batchSVMstruct(isubj).Y.SVML1 = YL1;
            batchSVMstruct(isubj).T.SVML1 = TL1;
            
            batchSVMstruct(isubj).null.AUCglm = AUCglmnull;
            batchSVMstruct(isubj).null.Xglm = Xglmnull;
            batchSVMstruct(isubj).null.Yglm = Yglmnull;
            batchSVMstruct(isubj).null.Tglm = Tglmnull;
            
            batchSVMstruct(isubj).mdlLM = mdlLM;
            batchSVMstruct(isubj).null.mdlLM = mdlLMnull;
            batchSVMstruct(isubj).mdlLMmin = mdlLM;
            batchSVMstruct(isubj).null.mdlLMmin = mdlLMnull;
            batchSVMstruct(isubj).minNodeType = minNodeType;
        end
        end
    end
    
    
    
    for i=1:Nsubj
        if ~isempty(batchSVMstruct(i).AUC)
        AUCas(i,1) = batchSVMstruct(i).AUC.SVM(1);
        AUCas(i,2) = batchSVMstruct(i).AUC.SVMallNodes(1);
        AUCas(i,3) = batchSVMstruct(i).AUC.SVMns(1);
        AUCas(i,4) = batchSVMstruct(i).null.AUC(1);
        end
    end
    
    % for i=1:27
    %     cAUC = batchSVMstruct(i).AUC.SVMminNode;
    %     cAUCnull = batchSVMstruct(i).null.AUC;
    %     AUCas(i,1:3) = cAUC;
    %     AUCas(i,4:6) = cAUCnull;
    % end
    
    if ~cvpart
        figure;
        boxplot(AUCas)
        set(gcf, 'color', 'w')
        box off
        if ~minNodeType
            xticklabels({'Selected Node Features' 'All Node Features' 'Non-selected Node Features' 'Random RT Label Shuffle'});
        elseif minNodeType==1
            xticklabels({'Feature Nodes' 'All Nodes' 'NS Nodes' 'Null'});
        else
            xticklabels({'Feature Nodes' 'All Nodes' 'NS Nodes' 'Null'});
        end
        ylabel('AUC')
        title(['SVM performance, thresh = ' num2str(percThresh)]);
    else
        if length(thresh) == 1
            figure; hold on
            mAUCasCVpart_inf = squeeze(mean(AUCasCVpart_inf,1));
            mAUCasCVpart_min = squeeze(mean(AUCasCVpart_min,1));
            mAUCasCVpart_all = squeeze(mean(AUCasCVpart_all,1));
            mAUCasCVpart_null = squeeze(mean(AUCasCVpart_null,1));
            
            errorbar(cvpart,mAUCasCVpart_inf(:,1),mAUCasCVpart_inf(:,1)-mAUCasCVpart_inf(:,2),mAUCasCVpart_inf(:,3)-mAUCasCVpart_inf(:,1), 'linewidth',1,'marker','.','markersize',18);
            errorbar(cvpart,mAUCasCVpart_min(:,1),mAUCasCVpart_min(:,1)-mAUCasCVpart_min(:,2),mAUCasCVpart_min(:,3)-mAUCasCVpart_min(:,1),'linewidth',1,'marker','.','markersize',18);
            errorbar(cvpart,mAUCasCVpart_all(:,1),mAUCasCVpart_all(:,1)-mAUCasCVpart_all(:,2),mAUCasCVpart_all(:,3)-mAUCasCVpart_all(:,1),'linewidth',1,'marker','.','markersize',18);
            errorbar(cvpart,mAUCasCVpart_null(:,1),mAUCasCVpart_null(:,1)-mAUCasCVpart_null(:,2),mAUCasCVpart_null(:,3)-mAUCasCVpart_null(:,1),'linewidth',1,'marker','d','markersize',12);
            legend('Influential nodes', 'Min/Max Slope Nodes', 'All Nodes', 'Null Model')
            title('SVM Performance on Remaining Test Data Based on Fractional Training Partition')
            xlabel('Fractional Training Partition (first n% trials)')
            ylabel('Mean AUC across subjects +/- 95% CI')
            set(gca,'FontSize',14)
            set(gcf,'color','w')
        end 
    end
    
    
    % for i=1:27
    %     rSq(i,1) = batchSVMstruct(i).mdlLM.Rsquared.Ordinary;
    %     rSq(i,2) = batchSVMstruct(i).mdlLM.Rsquared.Adjusted;
    %     rSq(i,3) = batchSVMstruct(i).null.mdlLM.Rsquared.Ordinary;
    % %     if ~useHalf
    % %     rSq(i,3) = batchSVMstruct(i).mdlLM.Rsquared.Ordinary*sum(behavStruct(i).ifast|behavStruct(i).islow);
    % %     else
    % %         rSq(i,3) = batchSVMstruct(i).mdlLM.Rsquared.Ordinary*sum(behavStruct(i).vRT>0&behavStruct(i).vRT<999);
    % %
    % %     end
    %     rSq(i,4) = batchSVMstruct(i).null.mdlLM.Rsquared.Adjusted;
    %
    %     rSq(i,5) = batchSVMstruct(i).mdlLMmin.Rsquared.Ordinary;
    %     rSq(i,6) = batchSVMstruct(i).mdlLMmin.Rsquared.Adjusted;
    %     rSq(i,7) = batchSVMstruct(i).null.mdlLMmin.Rsquared.Ordinary;
    % %     if ~useHalf
    % %     rSq(i,3) = batchSVMstruct(i).mdlLM.Rsquared.Ordinary*sum(behavStruct(i).ifast|behavStruct(i).islow);
    % %     else
    % %         rSq(i,3) = batchSVMstruct(i).mdlLM.Rsquared.Ordinary*sum(behavStruct(i).vRT>0&behavStruct(i).vRT<999);
    % %
    % %     end
    %     rSq(i,8) = batchSVMstruct(i).null.mdlLMmin.Rsquared.Adjusted;
    %
    % end
    
    % figure;
    % boxplot(rSq)
    % set(gcf, 'color', 'w')
    % box off
    % xticklabels({'Behav Ord R2' 'Null Ord R2' 'Behav Adj R2' 'Null Adj R2' 'Behav Min Ord R2' 'Null Min Ord R2' 'Behav Min Adj R2' 'Null Min Adj R2'});
    % title('Linear Regression against Continuous RT')
    
%     for i=1:BTRChange
%         for j=1:4
%             AUCf(i,j) = SVMstruct(i).fbands(j).AUC(1);
%         end
%     end
%     m = mean(AUCf,1);
%     for j=1:4
%         CI95(j,:) = vb_95CI(AUCf(:,j));
%     end
%     figure
%     errorbar(1,m(:,1),CI95(1,1),CI95(1,2))
%     hold on
%     errorbar(2,m(:,2),CI95(2,1),CI95(2,2))
%     errorbar(3,m(:,3),CI95(3,1),CI95(3,2))
%     errorbar(4,m(:,4),CI95(4,1),CI95(4,2))
%     bar(1,m(:,1))
%     bar(2,m(:,2))
%     bar(3,m(:,3))
%     bar(4,m(:,4))
%     box off
%     title('SVM per freq')
%     xticks([1:4])
%     xticklabels({'Theta/Alpha', 'Beta', 'LowGamma', 'HighGamma'})
%     set(gca,'FontSize',14)
%     set(gcf, 'color', 'w')
%     
%     figure
%     boxplot(AUCf)
%     box off
%     title('SVM per freq')
%     xticks([1:4])
%     xticklabels({'Theta/Alpha', 'Beta', 'LowGamma', 'HighGamma'})
%     set(gca,'FontSize',14)
%     set(gcf, 'color', 'w')
        
    if saveon
        save([svdir svnm,'.mat'],'batchSVMstruct', 'AUCas*', 'minNodeType', 'total*', 'LT*')
    end
    % end
    
%     close all
% AUCas(i,1) = batchSVMstruct(i).AUC.SVM(1);
AUC80_thresh(:,iThresh) = AUCas(:,1);
AUC80_thresh_null(:,iThresh) = AUCas(:,4);

totalSigNodesQ_all{iThresh} = totalSigNodesQ_thresh;
totalSigNodesP_all{iThresh} = totalSigNodesP_thresh;
totalSigNodesQmv_all{iThresh} = totalSigNodesQmv_thresh;
totalSigNodesPmv_all{iThresh} = totalSigNodesPmv_thresh;
if nullType>0
    totalSigNodesQnull_all{iThresh} = totalSigNodesQnull_thresh;
    totalSigNodesPnull_all{iThresh} = totalSigNodesPnull_thresh;
end
% totalSigNodesQnull_all{iThresh} = totalSigNodesQnull_thresh;
end

if length(thresh) >1
    %threshold selection figures
    figure;
    set(gcf,'Color','white')
    errorbar(thresh(1:end-1),mean(AUC80_thresh(:,1:end-1)),1.96.*std(AUC80_thresh(:,1:end-1))./sqrt(Nsubj),'-o',"color",'b',"MarkerSize",8,"MarkerFaceColor",'b')
    xlim([-.05,1])
    ylabel('Mean AUC Across Subjects +/- 95% CI')
    xlabel('Feature Selection Threshold (% of Repetitions Feature was Selected)')
    title('SVM Performance with Bootstrapped Feature Selection wtih 80% of Trials')
    
    %null comparison threshold selection
    figure;
    set(gcf,'Color','white')
    errorbar(thresh(1:end-1),mean(AUC80_thresh_null(:,1:end-1)),1.96.*std(AUC80_thresh_null(:,1:end-1))./sqrt(Nsubj),'-o',"color",'b',"MarkerSize",8,"MarkerFaceColor",'b')
    xlim([-.05,1])
    ylabel('Mean AUC Across Subjects +/- 95% CI')
    xlabel('Feature Selection Threshold (% of Repetitions Feature was Selected)')
    title('Null Model SVM Performance with Bootstrapped Feature Selection wtih 80% of Trials')
end

% AUCs_intratrial_pretrialFeats= AUCas(:,1);
% save(['/mnt/sdb1/CCDT/CCDTscripts/Figures/NullComparisonAUCs/AUCs_intratrial_pretrialFeats' num2str(nodeType) '.mat'],'AUCs_intratrial_pretrialFeats')
if saveForPlot
    save(['/mnt/sdb1/CCDT/CCDTscripts/Figures/Threshold_Selection/' plotSvNm],'AUC80_thresh','thresh'); %AUC80_thresh_pretrial_10
end

end % nodeType = r
%%
% for j=1:20
%     LTnullpow = LTnullStructRT20(j).nullpow;
%     LTnullqexp = LTnullStructRT20(j).nullqexp;
%     for i=1:BTRChange
%         x(i,j) = length(find(LTnullpow{i}(:,:,4)<0.05));
%         y(i,j) = length(find(LTnullqexp{i}(:,:,4)<0.05));
%         z(i,j) = length(find(LTpow{i}(:,:,4)<0.05));
%     end
% end
% xx = sum(x,1);
% yy = sum(y,1);
% zz = sum(y,1);
% totalSum = xx+yy;
% figure;
% histfit(totalSum,10)
% % hold on
% % plot([1452 1452], [0 19], 'k--', 'linewidth',2)
% 
% 
% for i=1:27
% x(i) = length(find(LTpow{i}(:,:,4)<0.05));
% y(i) = length(find(LTnullpow{i}(:,:,4)<0.05));
% end
% sum(x)
% sum(y)