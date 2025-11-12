function CCDT_graph_features
% calculate trial-by-trial graph communicability for CCDT cohort
% Correlate graph metrics with RT performance
%   VB 03/2020
% parameters
odir = 'test_output\'; % output directory
svnm = 'graph_features'; %name of output file
iev = 1; % event (1 = trial start cue, 2 = go cue, 3 = response)
win = [-500 0]; % perievent window (ms)
shiftDat = 0; %1 = create randomly shifted null data, 0= use regular data
fbands = [3 12; 12 30; 30 55; 70 110]; % frequency band ranges (Hz)
p.ddir = 'sample_raw_data\'; %'/CCDT/eeg/'; % raw data directory
p.subj = []; % subject (leave empty to batch process all subjects in database)
p.sess = []; % session date (leave empty to choose automatically)
p.stime = []; % session time (leave empty for all session times)
stsubj = 1; % start subject
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
cbehav = behav_sIall; %get behavioral data from all task sessions

%initiate structures to save
LTqexp = cell(Nsubj,1);
NFstruct = struct;
behavStruct = struct; % is this needed?

for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    p.subj = db{isubj,1}; p.sess = db{isubj,2};
    % load data
    [dat,Nsamp,fs,gch, gchlbl] = loadCCDTdata(p);
    
    % check contact labels to make sure they align with patient_loc file
    if length(gchlbl) == length(patient_loc(1).session(isubj).names)
        disp("gchs # = # chs in patient loc file")
    else
        disp("# gchs NOT the same as patient loc file")
        for ichlb = 1:length(patient_loc(1).session(isubj).names)
            for ichar = 1: length(patient_loc(1).session(isubj).names{ichlb})
                if gchlbl{ichlb}(ichar)~=patient_loc(1).session(isubj).names{ichlb}(ichar)
                    disp(['char mismatch in ch ' num2str(ichlb) '. channel removed.'])
                    gchlbl(ichlb) = [];
                    gch(ichlb) = [];
                    dat(:,ichlb) = [];
                end
            end
        end
        if length(gchlbl) == length(patient_loc(1).session(isubj).names)
            disp("gchs # = # chs in patient loc file")
        end
    end

    %remove contacts not labeled as gm or wm (outside of brain)
    gchlbl(patient_loc(1).session(isubj).type==0)=[];
    gch(patient_loc(1).session(isubj).type==0)=[];
    dat(:,patient_loc(1).session(isubj).type==0) = [];
    Nch = size(dat,2);
    
    if shiftDat %randomly shift data in time series as null comparison
        ndat = zeros(size(dat,1),size(dat,2));
        for i=1:size(dat,2)
        ndat(:,i) = circshift(dat(:,i),randi(size(dat,1)));
        end
        dat = ndat;
    end

    % load task events
    ccdt = parseCCDTevents(p);
    if ~isempty(Nsamp) && iscell(ccdt)
        nccdt = ccdt{1};
        disp(['Concatenating within-day multisession'])
        for ii = 2:length(ccdt)
            cccdt = ccdt{ii};
            cccdt(:,1:3) = cccdt(:,1:3) + sum(Nsamp(1:ii-1))*ones(size(cccdt,1),3); % event index accounting for concatenated data across session times
            nccdt = [nccdt; cccdt];
        end
        ccdt = nccdt;
    end
    Ntrl = size(ccdt,1); % # trials
    RT = ccdt(:,5); % reaction time (s)
    
    % channel pairs
    cgch = 1:Nch; % good channels
    chprs = nchoosek(cgch,2); % all pairs
    Npr = size(chprs,1); % # pairs
    disp(['Npr = ' num2str(Npr)]);
    
    % window indices for each trial relative to event/cue index
    swin = round(win(1)/1000*fs):round(win(2)/1000*fs);
    Nsamp = length(swin); % # samples in window
    indmat = ccdt(:,iev)*ones(1,Nsamp) + ones(Ntrl,1)*swin; % matrix of indices for all trial windows

    for ii = 1:Nbands % calculate features for each frequency band
        % bandpass filter to fband
        filtOrder = round(4/mean(fbands(ii,:))*fs); % length = 4 cycles of center frequency
        filtPts = fir1(filtOrder,fbands(ii,:)/(fs/2));
        fdat = filter(filtPts,1,dat);
        
        % instantaneous phase - Hilbert transform
        fdat = angle(hilbert(fdat));
        
        % window data by trial
        datwin = zeros(Nch,Ntrl,Nsamp);
        for jj = 1:Nch
            cdat = fdat(:,jj);
            datwin(jj,:,:) = cdat(indmat);
        end
        
        % calculate functional connectivity via single-trial PLV
        stplv = zeros(Npr,Ntrl);
        for jj = 1:Npr % calculate for all channel pairs
            x = squeeze(datwin(chprs(jj,1),:,:)); % trials x samples
            y = squeeze(datwin(chprs(jj,2),:,:));
            stplv(jj,:) = abs(mean(exp(1i*(x-y)),2)); % expected value across samples in window
        end

        % reorganize pairwise stplv into functional connectivity matrix
        Wst = zeros(Nch,Nch,Ntrl); 
        for jj = 1:Ntrl
            for kk = 1:Npr
                Wst(chprs(kk,1),chprs(kk,2),jj) = stplv(kk,jj);
                Wst(chprs(kk,2),chprs(kk,1),jj) = stplv(kk,jj);
            end
        end

        %set diagonal to 0 from NaN
        for i=1:size(Wst,1)
            for j=1:size(Wst,3)
                Wst(i,i,j)=0; 
            end
        end

        % calculate communicability for each subject per trial
        for i=1:size(Wst,3)
            NFqexp(:,:,i) = getCommunicability(Wst(:,:,i),1,1); %calculate pairwise communicability (Qexp)
        end
        NFqexp_nodal = squeeze(mean(NFqexp,1)); % expected value for each node
        
        % flag trials with bad RTs 
        inr = (RT<0 | RT>999); % incorrect responses, logical index
        disp([num2str(sum(inr)) ' trials omitted from analysis']);
        cTr = (~inr); %remove incorrect trials
        
        % regression - find features with a significant correlation with RT
        % includes all trials in this regression
        for ix=1:length(gch)
            [bqexp(ix,:),statsqexp(ix,:)] = robustfit(NFqexp_nodal(ix,cTr),RT(cTr)); 
            LTqexp{isubj}(ix,ii,:) = [bqexp(ix,2) bqexp(ix,2)-1.96*statsqexp(ix).se(2) bqexp(ix,2)+1.96*statsqexp(ix).se(2) statsqexp(ix).p(2)];
        end
               
        NFstruct(isubj).fbands(ii).NFqexp_nodal = NFqexp_nodal;
        NFstruct(isubj).gch = gch;
        NFstruct(isubj).gchlbl = gchlbl;
        
        behavStruct(isubj).subj = db{isubj,1};
        behavStruct(isubj).sess = db{isubj,2};
        behavStruct(isubj).vRT = RT;
        behavStruct(isubj).ifast = cbehav(isubj).ifast; % trials with RTs in the fastest 1/3
        behavStruct(isubj).islow = cbehav(isubj).islow; % trials with RTs in the slowest 1/3
        behavStruct(isubj).error = inr;
  
        clear NFq* bqexp stats*
    end
end

if saveon
save([odir svnm,'.mat'],'LT*','NFstruct', 'behavStruct')
end
end