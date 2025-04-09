function CCDTgraphRT2
% function CCDTgraphRT
% Correlate graph metrics with RT performance
%   VB 03/2020

% parameters
pdir = '/mnt/sdb1/CCDT/procData072023/'; % processed data directory
regr = 1; % robust regression? (0 or 1)
sI = 1; %0 for sI0, 1 for sIall
svnm = 'graphRT_040425_CAR_hg70110'; %'graphRT_sIall_500mspre_preextractionNull2'
iev = 1 % event (1 = trial start cue, 2 = go cue, 3 = response)
win = [-500 0] % perievent window (ms)
shiftDat = 0; %1 = create randomly shifted null data, 0= use regular data
useMid = 0; %use only middle third of trials for linear regression
fbands = [3 12; 12 30; 30 55; 70 110];
p.ddir = '/mnt/sdb1/CCDT/eeg/'; % data directory
p.subj = []; % subject (leave empty to batch process all subjects in database)
p.sess = []; % session date (leave empty to chose automatically)
p.stime = []; % session time (leave empty for all session times)
p.rln = 1; % remove line noise? (0 or 1)
p.rrf = 1 % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
p.outl = 1; % remove outlier (disconnected) channels? (0 or 1)
p.outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
p.outlThresh = 5; % outlier threshold (standard deviations)
stsubj = 1; % start subject
usePLI = 0; %0 = plv; 1 = pli
usePLV = 0; %0 = qexp; 1 = PLV/PLI (AS THE GRAPH ANALYSIS)
useSS = 0; % spatial smoothing: 0 = none, 1 = gaussian, 2 = neighbors
sigPThresh = 0.05;
saveon=1;


% load database
db = CCDTdatabase;
if ~isempty(p.subj)
    ind = cellfun(@(x) strcmp(x,p.subj),db(:,1)); % find specified subj
    if ~isempty(p.sess)
        ind2 = cellfun(@(x) strcmp(x,p.sess),db(:,2)); % find specified session
        ind = ind & ind2;
    end
    ind = find(ind);
    if isempty(ind), error('subj/sess not found in database'); end
    db = db(ind(1),:);
end

% load data
Nsubj = size(db,1); Nbands = size(fbands,1);
load('/home/mhedlund/Documents/RNET-main/CCDTscripts_BTR/allSubjBehaviorStruct.mat','behav','behav_sIall')
load('/home/mhedlund/Documents/RNET-main/CCDTscripts_BTR/patient_loc_101421.mat')
if sI
    cbehav = behav_sIall;
else
    cbehav = behav;
end
load('/home/mhedlund/Documents/RNET-main/CCDTscripts_BTR/all_WM.mat')


LTmod = zeros(Nsubj,size(fbands,1),4); % learning statistic [m +/- 95% ci, p]
LTstr = zeros(Nsubj,size(fbands,1),4);
LTqexp = cell(Nsubj,1);
sP = cell(Nsubj,size(fbands,1)); %index value within each subject good channels vector for most significant to least significant nodal Qexp RT prediction
sPval = cell(Nsubj,size(fbands,1));
NFstruct = struct;
SVMstruct = struct;
behavStruct = struct;


% AcrossTrialPLV = cell(Nsubj,1);
for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    p.subj = db{isubj,1}; q.subj = p.subj;
    p.sess = db{isubj,2}; q.sess = p.sess;
    q.ddir = [p.ddir(1:end-4) 'events/'];
    % load data
    [dat,Nsamp,fs,gch, gchlbl] = loadCCDTdata(p);
    % check contact labels to make sure they align with patient_loc file
    if length(gchlbl) == length(patient_loc(p.rrf).session(isubj).names)
        disp("gchs # = # chs in patient loc file")
    else
        disp("!!!! # gchs NOT the same as patient loc file")
        for ichlb = 1:length(patient_loc(p.rrf).session(isubj).names)
            for ichar = 1: length(patient_loc(p.rrf).session(isubj).names{ichlb})
                if gchlbl{ichlb}(ichar)~=patient_loc(p.rrf).session(isubj).names{ichlb}(ichar)
                    disp(['char mismatch in ch ' num2str(ichlb) '. channel removed.'])
                    gchlbl(ichlb) = [];
                    gch(ichlb) = [];
                    dat(:,ichlb) = [];
                end
            end
        end
        if length(gchlbl) == length(patient_loc(p.rrf).session(isubj).names)
            disp("gchs # = # chs in patient loc file")
        end
    end

    gchlbl(patient_loc(p.rrf).session(isubj).type==0)=[];
    gch(patient_loc(p.rrf).session(isubj).type==0)=[];
    dat(:,patient_loc(p.rrf).session(isubj).type==0) = [];
    Nch = size(dat,2);
    
    if useSS==1
        dat = fgsmooth(dat,10);
        
        
    elseif useSS==2
        datSS = zeros(size(dat));
        for i=2:size(dat,2)-1
           if strcmp(gchlbl{i}(1:2),gchlbl{i-1}(1:2))&&strcmp(gchlbl{i}(1:2),gchlbl{i+1}(1:2)) %both neighbors
            datSS(:,i) = 0.8.*dat(:,i) + 0.1.*dat(:,i-1) + 0.1.*dat(:,i+1);
               
           elseif strcmp(gchlbl{i}(1:2),gchlbl{i-1}(1:2))&&~strcmp(gchlbl{i}(1:2),gchlbl{i+1}(1:2)) %lower neighbor only
               datSS(:,i) = 0.8.*dat(:,i) + 0.2.*dat(:,i-1);
               
               
           elseif ~strcmp(gchlbl{i}(1:2),gchlbl{i-1}(1:2))&&strcmp(gchlbl{i}(1:2),gchlbl{i+1}(1:2)) %upper neighbor only
               datSS(:,i) = 0.8.*dat(:,i) + 0.2.*dat(:,i+1);
               
           end
            
        end
        datSS(:,1) = 0.8.*dat(:,1) + 0.2.*dat(:,2);
        datSS(:,end) = 0.8.*dat(:,end) + 0.2.*dat(:,end-1);
        dat = datSS;
        
    end
    
    if shiftDat %randomly shift data
        ndat = zeros(size(dat,1),size(dat,2));
        for i=1:size(dat,2)
        ndat(:,i) = circshift(dat(:,i),randi(size(dat,1)));
        end
        
        dat = ndat;
    end
    
    
    % load events
    ccdt = parseCCDTevents(q);
    if ~isempty(Nsamp) && iscell(ccdt)
        nccdt = ccdt{1};
        if sI
            disp(['Concatenating within-day multisession'])
            for ii = 2:length(ccdt)
                cccdt = ccdt{ii};
                cccdt(:,1:3) = cccdt(:,1:3) + sum(Nsamp(1:ii-1))*ones(size(cccdt,1),3); % event index accounting for concatenated data across session times
                nccdt = [nccdt; cccdt];
            end
        end
        ccdt = nccdt;
    end
    Ntrl = size(ccdt,1);
    CT = ccdt(:,1)/fs; % cue time (s)
    DT = ccdt(:,4); % delay time (s)
    RT = ccdt(:,5); % reaction time (s)
    
%     figure
%     plot(RT(RT>0&RT<999), '.')
%     title(['RT: subj ' num2str(isubj)])
    
    % channel pairs
    cgch = 1:Nch; % number of good channels (from loadCCDTdata)
    chprs = nchoosek(cgch,2); % all pairs
    Npr = size(chprs,1);
    disp(['Npr = ' num2str(Npr)]);
    
    % window indices
    swin = round(win(1)/1000*fs):round(win(2)/1000*fs);
    Nsamp = length(swin);
    indmat = ccdt(:,iev)*ones(1,Nsamp) + ones(Ntrl,1)*swin;
    
 
    
    
    for ii = 1:Nbands
        filtOrder = round(4/mean(fbands(ii,:))*fs); % length = 4 cycles of center frequency
        filtPts = fir1(filtOrder,fbands(ii,:)/(fs/2));
        fdat = filter(filtPts,1,dat);
        
        % instantaneous phase
        fdat = angle(hilbert(fdat));
        
        % window data
        datwin = zeros(Nch,Ntrl,Nsamp);
        for jj = 1:Nch
            cdat = fdat(:,jj);
            datwin(jj,:,:) = cdat(indmat);
        end
        
        % single-trial PLV
        stplv = zeros(Npr,Ntrl);
        stpli = zeros(Npr,Ntrl);
        for jj = 1:Npr
            x = squeeze(datwin(chprs(jj,1),:,:)); % trials x samples
            y = squeeze(datwin(chprs(jj,2),:,:));
            stplv(jj,:) = abs(mean(exp(1i*(x-y)),2)); % expected value across samples -- dp = angle(exp(1i*(x-y))); stplv2 = abs(mean(exp(1i*dp),2)); % alternative computation that gives same thing
%             stpli(jj,:) = abs(mean(signIm(x-y),2)); % from Stam et. al 2007 PLI, to get rid of artificial synchrony from volume conduction effects of the same signal or noise
        end
        
        % grouped trial PLV (FvS)
        iF = cbehav(isubj).ifast;
        iS = cbehav(isubj).islow;
        iL = cbehav(isubj).ilate;
        iE = cbehav(isubj).iearly;
        Ncond = 4;
        cond = false(Ncond,length(iF));
        cond(1,:) = iF;
        cond(2,:) = iS;
        cond(3,:) = iL;
        cond(4,:) = iE;
        atplv = zeros(Npr,Nsamp,Ncond);
        for jj = 1:Npr
            for ll=1:Ncond
                ccond = cond(ll,:);
            x = squeeze(datwin(chprs(jj,1),ccond,:)); % trials x samples 
            y = squeeze(datwin(chprs(jj,2),ccond,:));
            atplv(jj,:,ll) = abs(mean(exp(1i*(x-y)),1)); % expected value across trials
            end
        end
         
%         figure
%         plot(mean(atplv(:,:,1),1))
%         hold on
%         plot(mean(atplv(:,:,2),1))
%         plot(mean(atplv(:,:,3),1))
%         shg
%         plot(mean(atplv(:,:,4),1))
%         title(['Fband ' num2str(ii) ' subj: ' num2str(isubj)])
%         xlabel('Time')
%         ylabel('Averaged Wnole-Network Synchrony Across Trials (PLV)')
%         legend('Fast', 'Slow', 'Late', 'Early')
        
        % AcrossTrialPLV{isubj} = atplv;
        
        if usePLI
            stplv = stpli;
        end
      
        Wst = zeros(Nch,Nch,Ntrl);
        nodeStr = zeros(Nch,Ntrl);
        for jj = 1:Ntrl
            for kk = 1:Npr
                Wst(chprs(kk,1),chprs(kk,2),jj) = stplv(kk,jj);
                Wst(chprs(kk,2),chprs(kk,1),jj) = stplv(kk,jj);
            end
            nodeStr(1:Nch,jj) = sum(Wst(:,:,jj)); % nodal strength of desired network
        end
        
        NFstr = mean(nodeStr,1);
%         zNFstr = zscore(NFstr);
%         znodeStr = zscore(nodeStr);
        

        
        % Graph metrics for each subject per trial
        for i=1:size(Wst,1)
            for j=1:size(Wst,3)
                Wst(i,i,j)=0; %set diagonal to 0 from NaN
            end
        end
        for i=1:size(Wst,3)
            if ~usePLV
            NFqexp(:,:,i) = getCommunicability(Wst(:,:,i),1,1);
            else
               NFqexp(:,:,i) = Wst(:,:,i); 
            end
            [Ci(:,i),NFmod(:,i)] = community_louvain(Wst(:,:,i));
        end
        NFqexp_nodal = squeeze(mean(NFqexp,1));
        
%         if usePLV==2
%            NFqexp_nodal = NFmod'; 
%         end
%         
        for i=1:size(NFqexp,3)
           [Ciq(:,i),NFmodq(:,i)] = community_louvain(NFqexp(:,:,i)); 
        end
        
        % flag trials with bad RTs
        inr = (RT<0 | RT>999); % incorrect responses, logical index
        disp([num2str(sum(inr)) ' trials omitted from analysis']);
        
        if~useMid
            cTr = (~inr);
        else
            cTr = (~inr&~iF&~iS);
        end
        
                
        if regr
            for ix=1:length(gch)
                [bqexp(ix,:),statsqexp(ix,:)] = robustfit(NFqexp_nodal(ix,cTr),RT(cTr));
            end
            [bmod,statsmod] = robustfit(NFmod(cTr),RT(cTr));
            [bstr,statsstr] = robustfit(NFstr(cTr),RT(cTr));
            LTmod(isubj,ii,:) = [bmod(2) bmod(2)-1.96*statsmod.se(2) bmod(2)+1.96*statsmod.se(2) statsmod.p(2)];
            LTstr(isubj,ii,:) = [bstr(2) bstr(2)-1.96*statsstr.se(2) bstr(2)+1.96*statsstr.se(2) statsstr.p(2)];
            for ix=1:length(gch)
                LTqexp{isubj}(ix,ii,:) = [bqexp(ix,2) bqexp(ix,2)-1.96*statsqexp(ix).se(2) bqexp(ix,2)+1.96*statsqexp(ix).se(2) statsqexp(ix).p(2)];
            end
            [sPval{isubj,ii},sP{isubj,ii}] = sort(LTqexp{isubj}(:,ii,4), 'ascend');
        end
        
        % Differential Network Modularity
        try
            [Cidiff,NFmoddiff] = community_louvain(mean(Wst(:,:,iF),3)-mean(Wst(:,:,iS),3));
        catch
            Cidiff = []; NFmoddiff = [];
        end
        
        try
            [Cidiffq,NFmoddiffq] = community_louvain(mean(NFqexp(:,:,iF),3)-mean(NFqexp(:,:,iS),3));
        catch
            Cidiffq = []; NFmoddiffq = [];
        end
        
        %SVM
        cN = cell(size(Wst,3),1);
        cN(iF) = {'Fast'};
        cN(iS) = {'Slow'};
        cN(~iF&~iS) = [];
        
        try
            mdlSVMcv = fitcsvm(NFqexp_nodal((LTqexp{isubj}(:,ii,4)<=sigPThresh),(iF|iS))',cN,'KFold',5);
            SVMnodesUsed = gch(LTqexp{isubj}(:,ii,4)<=sigPThresh);
        catch
            mdlSVMcv = fitcsvm(NFqexp_nodal(sP{isubj,ii}(1),(iF|iS))',cN,'KFold',5);
            SVMnodesUsed = [];
            disp(['Freq ' num2str(ii) ' No sig nodes'])
        end
        [score_lbl,score_svm] = kfoldPredict(mdlSVMcv);
        [X,Y,T,AUC] = perfcurve(cN,score_svm(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
        
        % null comparison
        nullcN = cN(randperm(length(cN)));
        try
            mdlSVMcvnull = fitcsvm(NFqexp_nodal((LTqexp{isubj}(:,ii,4)<=sigPThresh),(iF|iS))',nullcN,'KFold',5);
        catch
            mdlSVMcvnull = fitcsvm(NFqexp_nodal(sP{isubj,ii}(1),(iF|iS))',nullcN,'KFold',5);
            disp(['Freq ' num2str(ii) ' No sig Null nodes'])
        end
        [score_lblnull,score_svmnull] = kfoldPredict(mdlSVMcvnull);
        [Xnull,Ynull,Tnull,AUCnull] = perfcurve(nullcN,score_svmnull(:,1),'Fast', 'NBoot', 1000, 'boottype', 'cper');
        
        
        %         NFstruct(isubj).fbands(ii).NFqexp = NFqexp;
       
            NFstruct(isubj).fbands(ii).NFqexp_nodal = NFqexp_nodal;
            NFstruct(isubj).fbands(ii).NFstr = NFstr;
            NFstruct(isubj).fbands(ii).NFmod = NFmod;
            NFstruct(isubj).fbands(ii).Ci = Ci;
            NFstruct(isubj).fbands(ii).NFmoddiff = NFmoddiff;
            NFstruct(isubj).fbands(ii).Cidiff = Cidiff;
            NFstruct(isubj).fbands(ii).NFmodq = NFmodq;
            NFstruct(isubj).fbands(ii).Ciq = Ciq;
            NFstruct(isubj).fbands(ii).NFmoddiffq = NFmoddiffq;
            NFstruct(isubj).fbands(ii).Cidiffq = Cidiffq;
            NFstruct(isubj).gch = gch;
            NFstruct(isubj).gchlbl = gchlbl;
            
            SVMstruct(isubj).fbands(ii).mdl = mdlSVMcv;
            SVMstruct(isubj).fbands(ii).nodes = SVMnodesUsed;
            SVMstruct(isubj).fbands(ii).score_lbl = score_lbl;
            SVMstruct(isubj).fbands(ii).score_svm = score_svm;
            SVMstruct(isubj).classNames = cN;
            SVMstruct(isubj).RT = RT(iF|iS);
            SVMstruct(isubj).fbands(ii).AUC = AUC;
            SVMstruct(isubj).fbands(ii).X = X;
            SVMstruct(isubj).fbands(ii).Y = Y;
            SVMstruct(isubj).fbands(ii).T = T;
            
            SVMstruct(isubj).fbands(ii).null.mdl = mdlSVMcvnull;
            SVMstruct(isubj).fbands(ii).null.score_lbl = score_lblnull;
            SVMstruct(isubj).fbands(ii).null.score_svm = score_svmnull;
            SVMstruct(isubj).fbands(ii).null.classNames = nullcN;
            SVMstruct(isubj).fbands(ii).null.AUC = AUCnull;
            SVMstruct(isubj).fbands(ii).null.X = Xnull;
            SVMstruct(isubj).fbands(ii).null.Y = Ynull;
            SVMstruct(isubj).fbands(ii).null.T = Tnull;
            
            behavStruct(isubj).subj = db{isubj,1};
            behavStruct(isubj).sess = db{isubj,2};
            behavStruct(isubj).vRT = RT;
            behavStruct(isubj).vDT = DT;
            behavStruct(isubj).ifast = iF;
            behavStruct(isubj).islow = iS;
            behavStruct(isubj).iearly = iE;
            behavStruct(isubj).ilate = iL;
            behavStruct(isubj).error = inr;
  
        
        
 

        
        
        clear NFq* NFstr NFmod* Ci* bmod bqexp bstr stats* mdlSVMcv cN score* X Y T AUC Xnull Ynull Tnull AUCnull SVMnodesUsed
        
        
    end
end

if saveon
save([pdir svnm,'.mat'],'LT*','sP*','NFstruct', 'SVMstruct', 'behavStruct')
end
end