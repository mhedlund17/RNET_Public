% CCDT_feature_details

% save details about selected features
% fPv --> mean z-scored value of selected power features over fastest third of trials
% fIDp --> fband of each selected power feature
% ssIDp --> subject/session number of selected power feature
% sPv --> mean z-scored value of selected power features over slowest third of trials
% aPv --> mean z-scored value of selected power features over all trials
% fQv, fIDq, ssIDq, sQv, aQv --> same variable meanings, but for qexp features

% parameters
% percent threshold: features with a significant regression slope in
% more than this percentage of bootstrap iterations will be included in
% the selected feature space. The threshold used for analysis in paper was
% chosen by evaluating SVM performance in predicting fast vs slow trials
% for 20 thresholds between 0 and 100%. 
percThresh = .3; % percentage, as a decimal
savedir = ''; %output directory
savefname = ''; %output file name
db = CCDTdatabase;
fbands = [3 12; 12 30; 30 55; 70 110]; 
Nsubj = size(db,1); Nbands = size(fbands,1)+1; Nfreq = height(fbands);

% load feature selection stats
fselect_dir = ''; % directory with feature selection file
fselect_fname = ''; % feature selection file name
load([fselect_dir fselect_fname],"sigPow","sigCom");

% load NFstruct for feature data and behavStruct for RT info
pdir = ''; % feature directory
pfname = ''; % power feature file name 
gfname = ''; % graph communicability feature file name
load([pdir pfname],"NFstruct"); 
NFstructP = NFstruct;
load([pdir gfname],"NFstruct","behavStruct");
NFstructQ = NFstruct;
clear NFstruct;

% identify features above selection threshold 
for isubj = 1:Nsubj
    sigNodesPow{isubj,1} = (sum(sigPow{isubj,1},3)/1000 > percThresh);
    sigNodesCom{isubj,1} = (sum(sigCom{isubj,1},3)/1000 > percThresh);
end

% initiate empty variables
fPv = []; sPv = []; %average power feature value over fast and slow  trials
fIDp = []; ssIDp = []; %frequency band and subjectID of power features
fQv = []; sQv = []; %fast and slow communicability features
aPv = []; aQv = []; %all power and communicability features
fIDq = []; ssIDq = []; %frequency band of communicability features
gchP = []; gchQ = []; %channel labels of power and communicability features
for isubj = 1:Nsubj
    disp(['Subject: ' num2str(isubj)])

    % get iF, iS - trial indices for fastest and slowest third of trials (by RT)
    RT = behavStruct(isubj).vRT;
    iF = behavStruct(isubj).ifast;
    iS = behavStruct(isubj).islow;

    %get features for current subject
    sigFeatP = sigNodesPow{isubj}; %nElec x nFreq
    sigFeatQ = sigNodesCom{isubj};
    nTrl = size(nPowDat,2);
    nFreq = size(nPowDat,3);

    for j = 1:nFreq %fbands, reformat features for current subject
        nComDat(:,:,j) = NFstructQ(isubj).fbands(j).NFqexp_nodal;
        nPowDat(:,:,j) = NFstructP(isubj).fbands(j).NFpow;
    end

    %get selected features
    for i = 1:nFreq
        %power features
        powFeaturesZ = zscore(nPowDat(:,:,i),1,'all'); % zscore over all channels/trials
        sigPowFeaturesZ = powFeaturesZ(sigFeatP(:,i),:); % get just selected channels
        fPv = vertcat(fPv,mean(sigPowFeaturesZ(:,iF),2)); %average over fast trials
        sPv = vertcat(sPv,mean(sigPowFeaturesZ(:,iS),2)); %average over slow trials
        aPv = vertcat(aPv,mean(sigPowFeaturesZ(:,:),2)); %average over all trials
        if ~isempty(sigPowFeaturesZ)
            fIDp = vertcat(fIDp,i*ones(height(sigPowFeaturesZ),1)); % frequency band
            ssIDp = vertcat(ssIDp,isubj*ones(height(sigPowFeaturesZ),1)); % subject ID
            gchP = vertcat(gchP,NFstructP(isubj).gchlbl(sigFeatP(:,i))); % channel label
        end
        %communicability features
        comFeaturesZ = zscore(nComDat(:,:,i),1,'all'); % zscore over all channels/trials
        sigComFeaturesZ = comFeaturesZ(sigFeatQ(:,i),:); % get just selected channels
        fQv = vertcat(fQv,mean(sigComFeaturesZ(:,iF),2)); %average over fast trials
        sQv = vertcat(sQv,mean(sigComFeaturesZ(:,iS),2)); %average over slow trials
        aQv = vertcat(aQv,mean(sigComFeaturesZ(:,:),2)); %average over all trials
        if ~isempty(sigComFeaturesZ)
            fIDq = vertcat(fIDq,i*ones(height(sigComFeaturesZ),1)); % frequency band
            ssIDq = vertcat(ssIDq,isubj*ones(height(sigComFeaturesZ),1)); % subject ID
            gchQ = vertcat(gchQ,NFstructQ(isubj).gchlbl(sigFeatQ(:,i))); % channel label
        end
    end
    clear nComDat nPowDat
end

% save details about selected features

save([savedir savefname],'fPv','sPv','fIDp','ssIDp','gchP','fQv','sQv','aQv','aPv','fIDq','ssIDq','gchQ')
