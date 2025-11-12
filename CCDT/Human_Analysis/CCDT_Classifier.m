function CCDT_Classifier(thresh)
% SVM classifier to predict if a trial is fast or slow
%   VB 03/2020
%   MH 08/2023 
warning off

%load features & basic info
pdir = 'test_output\'; % feature directory
pfname = 'power_features.mat'; % power feature file name
gfname = 'graph_features.mat'; % graph communicability feature file name
db = CCDTdatabase;
Nsubj = height(db);
fbands = [3 12; 12 30; 30 55; 70 110]; 

% load feature selection info
sfeat_dir = 'test_output\'; %feature selection directory
sfeat_fname = 'feature_selection.mat'; %feature selection file name - files saved in CCDTanalyze section 1

% save option, plot option, parameters
doROCPlot=0; %plots ROC
kfold=5; %for cross validation
percThresh = thresh/100; %percentage, as a decimal
saveon = 1; %save output
svdir = 'test_output\';
svnm = ['classifier_output_thresh' num2str(thresh)];

% load processed data
load([pdir pfname],"NFstruct"); 
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","behavStruct");
NFstructQ = NFstruct;
clear NFstruct;

% load selected features data 
load([sfeat_dir sfeat_fname],'sigPow','sigCom')  
for isubj = 1:Nsubj
    % identify features above selection threshold
    sigNodesPow{isubj,1} = (sum(sigPow{isubj,1},3)/1000 > percThresh);
    sigNodesCom{isubj,1} = (sum(sigCom{isubj,1},3)/1000 > percThresh);
end
clear sigPow* sigCom*
    
for isubj = 1:Nsubj        
    disp(['Subject: ' num2str(isubj)])

    %initiate empty structures
    qexpNodes = cell(size(fbands,1),3); %raw channel, channel index, channel label
    powNodes = cell(size(fbands,1),3);

    % get info for current subject
    cbehavStruct = behavStruct(isubj);
    cNFstructP = NFstructP(isubj);
    cNFstructQ = NFstructQ(isubj);
    cCHpow = cNFstructP.gch;
    cCH = cNFstructQ.gch;
    cgchlblpow = cNFstructP.gchlbl;
    cgchlbl = cNFstructQ.gchlbl;

    % setup behavioral labels
    cRT = cbehavStruct.vRT; % reaction time (ms)
    goodTrl = ~(cRT<50 | cRT>999); %remove incorrect trials
    iF = cbehavStruct.ifast; %fastest third
    iS = cbehavStruct.islow; %slowest third
    cN = cell(length(cRT),1);
    cN(iF) = {'Fast'};
    cN(iS) = {'Slow'};
    cN(~iF&~iS) = []; % remove labels except fastest/slowest third
    cNRT = cRT(iF|iS); % remove RTs except fastest/slowest third
               
    % null comparison: shuffled label
    ncN = cN(randperm(length(cN)));
           
    %frequency-specific features
    tNFpow = cNFstructP.fbands(1).NFpow; %theta/alpha band
    bNFpow = cNFstructP.fbands(2).NFpow; %beta band
    lgNFpow = cNFstructP.fbands(3).NFpow; %low gamma band
    hgNFpow = cNFstructP.fbands(4).NFpow; %high gamma band
    
    tNFqexp = cNFstructQ.fbands(1).NFqexp_nodal;
    bNFqexp = cNFstructQ.fbands(2).NFqexp_nodal;
    lgNFqexp = cNFstructQ.fbands(3).NFqexp_nodal;
    hgNFqexp = cNFstructQ.fbands(4).NFqexp_nodal;
    
    %get index for selected nodes
    clear useQnodes usePnodes 
    for i=1:length(fbands)
        useQnodes{i} = find(sigNodesCom{isubj}(:,i));
        usePnodes{i} = find(sigNodesPow{isubj}(:,i));
        nsQnodes{i} = find(~sigNodesCom{isubj}(:,i));
        nsPnodes{i} = find(~sigNodesPow{isubj}(:,i));
    end
    totalSigNodesQ_thresh(isubj,:) = cellfun(@length,useQnodes);
    totalSigNodesP_thresh(isubj,:) = cellfun(@length,usePnodes);
        
    % combine features at just fast/slow trials in all fbands for selected nodes, non-selected nodes, and all nodes
    dat = [tNFqexp(useQnodes{1},(iF|iS))' bNFqexp(useQnodes{2},(iF|iS))' lgNFqexp(useQnodes{3},(iF|iS))' hgNFqexp(useQnodes{4},(iF|iS))']; %selected features
    nsDat = [tNFqexp(nsQnodes{1},(iF|iS))' bNFqexp(nsQnodes{2},(iF|iS))' lgNFqexp(nsQnodes{3},(iF|iS))' hgNFqexp(nsQnodes{4},(iF|iS))']; %non-selected features
    datAllNodes = [tNFpow(:,(iF|iS))' bNFpow(:,(iF|iS))' lgNFpow(:,(iF|iS))' hgNFpow(:,(iF|iS))' tNFqexp(:,(iF|iS))' bNFqexp(:,(iF|iS))' lgNFqexp(:,(iF|iS))' hgNFqexp(:,(iF|iS))']; %all features
       
    %create SVM model for selected and non-selected features
    try
        mdlSVMcv = fitcsvm(dat,cN,'KFold',kfold);
        mdlSVMcvNS = fitcsvm(nsDat(:,randi(size(nsDat,2),size(dat,2))),cN,'KFold',kfold);
        mdlSVMcv_allNodes = fitclinear(datAllNodes,cN,'learner', 'svm','KFold',kfold,'regularization','ridge');
    catch
        disp(['No significant feature nodes for subject session ' num2str(isubj) ' -- skipping...'])
        break
    end

    % model prediction for selected features, non-selected features, all nodes
    [score_lbl,score_svm] = kfoldPredict(mdlSVMcv);
    [X,Y,T,AUC] = perfcurve(cN,score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');

    [score_lblNS,score_svmNS] = kfoldPredict(mdlSVMcvNS);
    [Xns,Yns,Tns,AUCns] = perfcurve(cN,score_svmNS(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');

    [~,score_svmAllNodes] = kfoldPredict(mdlSVMcv_allNodes);
    [Xall,Yall,Tall,AUCall] = perfcurve(cN,score_svmAllNodes(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
     
    % shuffled label null comparison model
    mdlSVMcvnull = fitcsvm(dat,ncN,'KFold',kfold);
    [score_lblnull,score_svmnull] = kfoldPredict(mdlSVMcvnull);
    [Xnull,Ynull,Tnull,AUCnull] = perfcurve(ncN,score_svmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
                         
    if exist("AUC")
        %store model info for selected features
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
        
        % store model info for shuffled label null
        batchSVMstruct(isubj).null.mdl = mdlSVMcvnull;
        batchSVMstruct(isubj).null.score_lbl = score_lblnull;
        batchSVMstruct(isubj).null.score_svm = score_svmnull;
        batchSVMstruct(isubj).null.classNames = ncN;
        batchSVMstruct(isubj).null.AUC = AUCnull;
        batchSVMstruct(isubj).null.X = Xnull;
        batchSVMstruct(isubj).null.Y = Ynull;
        batchSVMstruct(isubj).null.T = Tnull;
        
        % store model output for all nodes
        batchSVMstruct(isubj).AUC.SVMallNodes = AUCall;
        batchSVMstruct(isubj).X.SVMallNodes = Xall;
        batchSVMstruct(isubj).Y.SVMallNodes = Yall;
        batchSVMstruct(isubj).T.SVMallNodes = Tall;
        
        %store model output for non-selected nodes
        batchSVMstruct(isubj).AUC.SVMns = AUCns;
        batchSVMstruct(isubj).X.SVMns = Xns;
        batchSVMstruct(isubj).Y.SVMns = Yns;
        batchSVMstruct(isubj).T.SVMns = Tns;

        if doROCPlot
            % ROC plot single subject performance selected features vs. null
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
            title('Behavioral (blue) vs Shuffled Label Null (gray) ROC')
                
            %add ROC subplot for subject's performance selected features vs. null
            figure(2)
            set(gcf,'Color','white')
            sgtitle('Behavioral (blue) vs Shuffled Label Null (gray) ROC')
            subplot(3,9,isubj)
            hold on
            shadedErrorBar(X(:,1),Y(:,1),Y(:,1)-Y(:,2),'b');
            shadedErrorBar(Xnull(:,1),Ynull(:,1),Ynull(:,1)-Ynull(:,2),'k');
            title(['Subject ' num2str(isubj)])
            plot([0,1],[0,1],'r')
        end
    end
end %subj loop
    
% compile AUC's across subjects
for i=1:Nsubj
    if ~isempty(batchSVMstruct(i).AUC)
        AUCas(i,1) = batchSVMstruct(i).AUC.SVM(1);
        AUCas(i,2) = batchSVMstruct(i).AUC.SVMallNodes(1);
        AUCas(i,3) = batchSVMstruct(i).AUC.SVMns(1);
        AUCas(i,4) = batchSVMstruct(i).null.AUC(1);
    end
end
    
%plot AUC distributions across subjects
figure;
boxplot(AUCas)
set(gcf, 'color', 'w')
box off
xticklabels({'Feature Nodes' 'All Nodes' 'NS Nodes' 'Null'});
ylabel('AUC')
title(['SVM performance, thresh = ' num2str(percThresh)]);
          
if saveon
    save([svdir svnm,'.mat'],'batchSVMstruct', 'AUCas*', 'total*')
end
