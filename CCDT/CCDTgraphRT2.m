function CCDTgraphRT2
% function CCDTgraphRT
% calculate graph communicability
% Correlate graph metrics with RT performance
%   VB 03/2020

% parameters
odir = ''; % output directory
regr = 1; % robust regression? (0 or 1)
sI = 1; %0 for first task session only, 1 to use all task sessions in same day
svnm = ''; %name of output file
iev = 1 % event (1 = trial start cue, 2 = go cue, 3 = response)
win = [-500 0] % perievent window (ms)
shiftDat = 0; %1 = create randomly shifted null data, 0= use regular data
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
usePLV = 0; %0 = qexp; 1 = PLV node strength (AS THE GRAPH ANALYSIS)
useSS = 0; % spatial smoothing: 0 = none, 1 = gaussian, 2 = neighbors
saveon=1; %save features

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
load('allSubjBehaviorStruct.mat','behav','behav_sIall') %behavioral data
load('patient_loc_101421.mat') %anatomical labels
if sI
    cbehav = behav_sIall;
else
    cbehav = behav;
end
load('all_WM.mat')

LTqexp = cell(Nsubj,1);
sP = cell(Nsubj,size(fbands,1)); %index value within each subject good channels vector for most significant to least significant nodal Qexp RT prediction
sPval = cell(Nsubj,size(fbands,1));
NFstruct = struct;
SVMstruct = struct;
behavStruct = struct;

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
        disp("# gchs NOT the same as patient loc file")
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
    %remove contacts not labeled as gm or wm
    gchlbl(patient_loc(p.rrf).session(isubj).type==0)=[];
    gch(patient_loc(p.rrf).session(isubj).type==0)=[];
    dat(:,patient_loc(p.rrf).session(isubj).type==0) = [];
    Nch = size(dat,2);
    
    if useSS==1 %gaussian spatial smoothing
        dat = fgsmooth(dat,10);        
    elseif useSS==2 %neighbors spatial smoothing
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
    
    % channel pairs
    cgch = 1:Nch; % number of good channels
    chprs = nchoosek(cgch,2); % all pairs
    Npr = size(chprs,1);
    disp(['Npr = ' num2str(Npr)]);
    
    % window indices
    swin = round(win(1)/1000*fs):round(win(2)/1000*fs);
    Nsamp = length(swin);
    indmat = ccdt(:,iev)*ones(1,Nsamp) + ones(Ntrl,1)*swin;

    for ii = 1:Nbands
        %filter to fband
        filtOrder = round(4/mean(fbands(ii,:))*fs); % length = 4 cycles of center frequency
        filtPts = fir1(filtOrder,fbands(ii,:)/(fs/2));
        fdat = filter(filtPts,1,dat);
        
        % instantaneous phase
        fdat = angle(hilbert(fdat));
        
        % window data by trial
        datwin = zeros(Nch,Ntrl,Nsamp);
        for jj = 1:Nch
            cdat = fdat(:,jj);
            datwin(jj,:,:) = cdat(indmat);
        end
        
        % single-trial PLV
        stplv = zeros(Npr,Ntrl);
        for jj = 1:Npr
            x = squeeze(datwin(chprs(jj,1),:,:)); % trials x samples
            y = squeeze(datwin(chprs(jj,2),:,:));
            stplv(jj,:) = abs(mean(exp(1i*(x-y)),2)); % expected value across samples -- dp = angle(exp(1i*(x-y))); stplv2 = abs(mean(exp(1i*dp),2)); % alternative computation that gives same thing
        end
        
        % grouped trial PLV (Fast vs Slow)
        iF = cbehav(isubj).ifast;
        iS = cbehav(isubj).islow;
        Ncond = 2;
        cond = false(Ncond,length(iF));
        cond(1,:) = iF;
        cond(2,:) = iS;

        % calculate functional connectivity
        Wst = zeros(Nch,Nch,Ntrl); %functional connectivity matrix
        for jj = 1:Ntrl
            for kk = 1:Npr
                Wst(chprs(kk,1),chprs(kk,2),jj) = stplv(kk,jj);
                Wst(chprs(kk,2),chprs(kk,1),jj) = stplv(kk,jj);
            end
        end

        % Graph metrics for each subject per trial
        for i=1:size(Wst,1)
            for j=1:size(Wst,3)
                Wst(i,i,j)=0; %set diagonal to 0 from NaN
            end
        end
        for i=1:size(Wst,3)
            if ~usePLV %use nodal communicability as graph metric
                NFqexp(:,:,i) = getCommunicability(Wst(:,:,i),1,1);
            else %use plv nodal strength as graph metric
                NFqexp(:,:,i) = Wst(:,:,i); 
            end
        end
        NFqexp_nodal = squeeze(mean(NFqexp,1));
        
        % flag trials with bad RTs
        inr = (RT<0 | RT>999); % incorrect responses, logical index
        disp([num2str(sum(inr)) ' trials omitted from analysis']);
        cTr = (~inr); %remove incorrect trials
        
        % regression   
        for ix=1:length(gch)
            [bqexp(ix,:),statsqexp(ix,:)] = robustfit(NFqexp_nodal(ix,cTr),RT(cTr));
        end
        for ix=1:length(gch)
            LTqexp{isubj}(ix,ii,:) = [bqexp(ix,2) bqexp(ix,2)-1.96*statsqexp(ix).se(2) bqexp(ix,2)+1.96*statsqexp(ix).se(2) statsqexp(ix).p(2)];
        end
        [sPval{isubj,ii},sP{isubj,ii}] = sort(LTqexp{isubj}(:,ii,4), 'ascend');
               
        NFstruct(isubj).fbands(ii).NFqexp_nodal = NFqexp_nodal;
        NFstruct(isubj).gch = gch;
        NFstruct(isubj).gchlbl = gchlbl;
        
        behavStruct(isubj).subj = db{isubj,1};
        behavStruct(isubj).sess = db{isubj,2};
        behavStruct(isubj).vRT = RT;
        behavStruct(isubj).ifast = iF;
        behavStruct(isubj).islow = iS;
        behavStruct(isubj).error = inr;
  
        clear NFq* bqexp stats*
    end
end

if saveon
save([odir svnm,'.mat'],'LT*','sP*','NFstruct', 'behavStruct')
end
end