% function CCDTBatchClassifier2
% Train SVM classifier on combined neural data -- WM vs. GM with full
% network
%   VB 03/2020
%   MH 08/2023 
clear; close all; clc

warning off
% parameters
%load features
pdir = ''; % feature directory
pfname = ''; % power feature file name
gfname = ''; % graph communicability feature file name

%load behavioral data path
rtdir = '/mnt/sdb1/CCDT/orig_procData/'; % behavioral data directory
rtfnm = 'graphRT_sIall_500mspre'; %use graph metric versus PLV metric defined here

% load selected features
sdir = ''; %feature selection directory
fname8020 = ''; %feature selection file name

plotSvNm = 'AUC80_thresh_CAR_hg70110_intra_091124.mat';
saveForPlot = 0;
Qpolarity = 1; %1 = CAR, 2 = BP
Ppolarity = 1; %1 = CAR, 2 = BP
fbands = [3 12; 12 30; 30 55; 70 110]; 
saveon = 0;
kfold=5;
nullType = 0; % 0 for shuffle labels, 1 for simulated vRT, 2 for permuted vRT
usePrevTrial = 0; % if nullType == 0, use previous trial (1) or randomly shuffle (0) for post extraction null
cvpart = 0;%[0.05:0.05:0.85]; %0 for use full set with kfold, otherwise specify fraction for temporal partition
sigPThresh = 0.05; % P threshold for feature selection
db = CCDTdatabase;
% thresh = (0:5:100)/100; %set range of thresholds
thresh = .3; %set one threshold after selecting from range

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

load([pdir pfname],"NFstruct","LTpow"); %
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","LTqexp");
NFstructQ = NFstruct;
clear NFstruct;
load([rtdir rtfname],"behavStruct");%load behavStruct RT data
Nsubj = height(db);
load('all_WM.mat')
load patient_loc_120623 %anatomical labels

for iThresh = 1:length(thresh)
    percThresh = thresh(iThresh);
    disp(['Threshold: ' num2str(thresh(iThresh))])

    % load selected features data 
    load([dir8020 fname8020],'sigPow','sigCom')  
    for isubj = 1:Nsubj
        sigNodesPow{isubj,1} = (sum(sigPow{isubj,1},3)/1000 > percThresh);
        sigNodesCom{isubj,1} = (sum(sigCom{isubj,1},3)/1000 > percThresh);
    end
    if nullType == 2
        load([dir8020 fname8020],'nullShuffleRT','sigPowNRT','sigComNRT')
        for isubj = 1:Nsubj
            sigNodesPowNRT{isubj,1} = (sum(sigPowNRT{isubj,1},3)/1000 > percThresh);
            sigNodesComNRT{isubj,1} = (sum(sigComNRT{isubj,1},3)/1000 > percThresh);
        end 
    end
    clear sigPow* sigCom*

    LTnullStructRT20 = struct;
    
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

        % find correct trials
        cRT = cbehavStruct.vRT; % reaction time (ms)
        goodTrl = ~(cRT<50 | cRT>999);

        %Behavior setup
        iF = cbehavStruct.ifast;
        iS = cbehavStruct.islow;
        cN = cell(length(cRT),1);
        cN(iF) = {'Fast'};
        cN(iS) = {'Slow'};
        cN(~iF&~iS) = [];
        cNRT = cRT(iF|iS);
               
        if nullType==1
            nullRT = rand(length(cRT),1).*999;
            inr = (cRT<50 | cRT>999);
        elseif nullType==2
            nullRT = nullShuffleRT{isubj};
            inr = (nullRT<50 | nullRT>999);
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
            iFn(cB(1:round(sum(ncTc)/3))) = 1;
            iSn(cB(end-round(sum(ncTc)/3)+1:end)) = 1;
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
         
        [~,tPmSi] = min(LTpow{isubj}(:,1,1)); [~,tPmaSi] = max(LTpow{isubj}(:,1,1));
        [~,bPmSi] = min(LTpow{isubj}(:,2,1)); [~,bPmaSi] = max(LTpow{isubj}(:,2,1));
        [~,lgPmSi] = min(LTpow{isubj}(:,3,1)); [~,lgPmaSi] = max(LTpow{isubj}(:,3,1));
        [~,hgPmSi] = min(LTpow{isubj}(:,4,1)); [~,hgPmaSi] = max(LTpow{isubj}(:,4,1));
        
        [~,tQmSi] = min(LTqexp{isubj}(:,1,1)); [~,tQmaSi] = max(LTqexp{isubj}(:,1,1));
        [~,bQmSi] = min(LTqexp{isubj}(:,2,1)); [~,bQmaSi] = max(LTqexp{isubj}(:,2,1));
        [~,lgQmSi] = min(LTqexp{isubj}(:,3,1)); [~,lgQmaSi] = max(LTqexp{isubj}(:,3,1));
        [~,hgQmSi] = min(LTqexp{isubj}(:,4,1)); [~,hgQmaSi] = max(LTqexp{isubj}(:,4,1));
        
        clear useQnodes usePnodes 
        for i=1:length(fbands)
            useQnodes{i} = find(sigNodesCom{isubj}(:,i));
            usePnodes{i} = find(sigNodesPow{isubj}(:,i));
            nsQnodes{i} = find(~sigNodesCom{isubj}(:,i));
            nsPnodes{i} = find(~sigNodesPow{isubj}(:,i));
             
            if nullType==2
                nullNodesQ{i} = find(sigNodesComNRT{isubj}(:,i));
                nullNodesP{i} = find(sigNodesPowNRT{isubj}(:,i));
            end
        end
            
        totalSigNodesQ_thresh(isubj,:) = cellfun(@length,useQnodes);
        totalSigNodesP_thresh(isubj,:) = cellfun(@length,usePnodes);
        if nullType>0
            totalSigNodesQnull_thresh(isubj,:) = cellfun(@length,nullNodesQ);
            totalSigNodesPnull_thresh(isubj,:) = cellfun(@length,nullNodesP);
        end
        
        dat = [tNFqexp(useQnodes{1},(iF|iS))' bNFqexp(useQnodes{2},(iF|iS))' lgNFqexp(useQnodes{3},(iF|iS))' hgNFqexp(useQnodes{4},(iF|iS))']; %selected features
        nsDat = [tNFqexp(nsQnodes{1},(iF|iS))' bNFqexp(nsQnodes{2},(iF|iS))' lgNFqexp(nsQnodes{3},(iF|iS))' hgNFqexp(nsQnodes{4},(iF|iS))']; %non-selected features
        datAllNodes = [tNFpow(:,(iF|iS))' bNFpow(:,(iF|iS))' lgNFpow(:,(iF|iS))' hgNFpow(:,(iF|iS))' tNFqexp(:,(iF|iS))' bNFqexp(:,(iF|iS))' lgNFqexp(:,(iF|iS))' hgNFqexp(:,(iF|iS))']; %all features
        
        if nullType>0
            nullDat = [tNFqexp(nullNodesQ{1},(iFn|iSn))' bNFqexp(nullNodesQ{2},(iFn|iSn))' lgNFqexp(nullNodesQ{3},(iFn|iSn))' hgNFqexp(nullNodesQ{4},(iFn|iSn))'];
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
                catch
                    disp(['No significant feature nodes for subject session ' num2str(isubj) ' -- skipping...'])
                    break
                end
                [score_lbl,score_svm] = kfoldPredict(mdlSVMcv);
                [X,Y,T,AUC] = perfcurve(cN,score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                [score_lblNS,score_svmNS] = kfoldPredict(mdlSVMcvNS);
                [Xns,Yns,Tns,AUCns] = perfcurve(cN,score_svmNS(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                %logistic regression comparison
                mdlGLMcv = fitclinear(dat,cN,'learner','logistic','kfold',kfold,'regularization','ridge');
                [~,score_glm] = kfoldPredict(mdlGLMcv);
                [Xglm,Yglm,Tglm,AUCglm] = perfcurve(cN,score_glm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                mdlGLMcvNS = fitclinear(nsDat,cN,'learner','logistic','kfold',kfold,'regularization','ridge');
                [~,score_glmNS] = kfoldPredict(mdlGLMcvNS);
                [XglmNS,YglmNS,TglmNS,AUCglmNS] = perfcurve(cN,score_glmNS(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                % null comparison
                if ~nullType
                    mdlSVMcvnull = fitcsvm(dat,ncN,'KFold',kfold);
                    mdlGLMcvnull = fitclinear(dat,ncN,'learner','logistic','kfold',kfold, 'regularization','ridge');
                else
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
                
                % SVM with ridge and lasso
                mdlSVML2 = fitclinear(dat,cN,'learner','svm','regularization','ridge', 'kfold',kfold);
                [score_lblL2,score_svmL2] = kfoldPredict(mdlSVML2);
                [XL2,YL2,TL2,AUCL2] = perfcurve(cN,score_svmL2(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                mdlSVML1 = fitclinear(dat,cN,'learner','svm','regularization','lasso', 'kfold',kfold);
                [score_lblL2,score_svmL1] = kfoldPredict(mdlSVML1);
                [XL1,YL1,TL1,AUCL1] = perfcurve(cN,score_svmL1(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                
                % test all Nodes
                mdlSVMcv_allNodes = fitclinear(datAllNodes,cN,'learner', 'svm','KFold',kfold,'regularization','ridge');
                [~,score_svmAllNodes] = kfoldPredict(mdlSVMcv_allNodes);
                [Xall,Yall,Tall,AUCall] = perfcurve(cN,score_svmAllNodes(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
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
            end
        end
        
        if exist("AUC")
        AUCasCVpart_inf(isubj,:,:) = AUC;
        AUCasCVpart_all(isubj,:,:) = AUCall;
        AUCasCVpart_null(isubj,:,:) = AUCnull;
        AUCasCVpart_NS(isubj,:,:) = AUCns;
        
        batchSVMstruct(isubj).mdl.SVM = mdlSVMcv;
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
        
        %plot performance
        if length(thresh)==1
            %figure 2D
            figure(1)
            set(gcf,'Color','white','Units','inches','Position',[2,2,5,3.5])
            hold on;
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
    
    if ~cvpart
        figure;
        boxplot(AUCas)
        set(gcf, 'color', 'w')
        box off
        xticklabels({'Feature Nodes' 'All Nodes' 'NS Nodes' 'Null'});
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
        
    if saveon
        save([svdir svnm,'.mat'],'batchSVMstruct', 'AUCas*', 'total*', 'LT*')
    end

    AUC80_thresh(:,iThresh) = AUCas(:,1);
    AUC80_thresh_null(:,iThresh) = AUCas(:,4);
    
    totalSigNodesQ_all{iThresh} = totalSigNodesQ_thresh;
    totalSigNodesP_all{iThresh} = totalSigNodesP_thresh;
    if nullType>0
        totalSigNodesQnull_all{iThresh} = totalSigNodesQnull_thresh;
        totalSigNodesPnull_all{iThresh} = totalSigNodesPnull_thresh;
    end
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

if saveForPlot
    save([svdir plotSvNm],'AUC80_thresh','thresh');
end