% CCDT_trPLV
% calculate time resolved phase locking value (PLV) for CCDT cohort
% Time-resolved PLV is calculated by averaging phase locking over a set of
% trials corresponding to the following conditions:
%   1. Fastest 1/3 of trials
%   2. Slowest 1/3 of trials
%   3. All non-error trials
%   4. Middle 1/3 of trials
%   MH 07/2025
clc; clear; close all 

% parameters
odir = ''; % output directory
svnm = ''; %name of output file
iev = 1; % event (1 = trial start cue, 2 = go cue, 3 = response)
win = [-2100 600]; % perievent window (ms) --> iev=1: [-1100, 1600]; iev=2: [-2100, 600];
fbands = [3 12; 12 30; 30 55; 70 110]; % frequency band ranges (Hz)
p.ddir = '/sample_raw_data/'; % raw data directory
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

% load anatomical and behavioral data
Nsubj = size(db,1); Nbands = size(fbands,1);
load('allSubjBehaviorStruct.mat','behav','behav_sIall')
load('patient_loc_101421.mat')
cbehav = behav_sIall;
load('/home/mhedlund/Documents/RNET-main/CCDTscripts_BTR/all_WM.mat')

%initiate structures to save
TRstruct = struct;
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

    %remove contacts not labeled as gm or wm (outside of brain)
    gchlbl(patient_loc(p.rrf).session(isubj).type==0)=[];
    gch(patient_loc(p.rrf).session(isubj).type==0)=[];
    dat(:,patient_loc(p.rrf).session(isubj).type==0) = [];
    Nch = size(dat,2);
    
    % load task events
    ccdt = parseCCDTevents(q);
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

    % flag trials with bad RTs
    inr = (RT<0 | RT>999); % incorrect responses, logical index    
    disp([num2str(sum(inr)) ' trials omitted from analysis']);
   
    for ii = 1:Nbands
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
 
        % grouped trial PLV (FvS)
        iF = cbehav(isubj).ifast; % fastest 1/3 of trials
        iS = cbehav(isubj).islow; % slowest 1/3 of trials
        condLabels = {'fast','slow','all','mid'}; % condition labels
        Ncond = size(condLabels,2); % # conditions
        cond = false(Ncond,length(iF));
        cond(1,:) = iF;
        cond(2,:) = iS;
        cond(3,:) = ~inr; %all trials, but remove invalid RTs (these are already removed from iF and iS)
        cond(4,:) = ~iF & ~iS & ~inr; % middle 1/3

        % calculate time-resolved plv across trials
        atplv = zeros(Npr,Nsamp,Ncond);
        for jj = 1:Npr %calculate for all channel pairs
            for ll=1:Ncond
                ccond = cond(ll,:);
                x = squeeze(datwin(chprs(jj,1),ccond,:)); % samples x trials
                y = squeeze(datwin(chprs(jj,2),ccond,:));
                cNtrl = size(x,2); %number of trials in this condition
                atplv(jj,:,ll) = abs(mean(exp(1i*(x-y)),1)); % expected value across trials, Npr x Nsamp x Ncond
            end
        end

        %save time-resolved data
        TRstruct.fbands(ii).atplv = atplv;
        TRstruct.Nch = Nch;
        TRstruct.Nsamp = Nsamp;
        TRstruct.Npr = Npr;
        TRstruct.iev = iev;
        TRstruct.win = win;
        TRstruct.winSamps = swin;
        TRstruct.t = timeaxis;
        TRstruct.fs = fs;
        TRstruct.gch = gch;
        TRstruct.gchlbl = gchlbl;
    end
    
    %save behavioral info
    behavStruct.subj = db{isubj,1};
    behavStruct.sess = db{isubj,2};
    behavStruct.vRT = RT;
    behavStruct.vDT = DT;
    behavStruct.ifast = iF;
    behavStruct.islow = iS;
    behavStruct.error = inr;
    behavStruct.cond = cond;
    behavStruct.condLabels = condLabels;
    
    if saveon
        save([odir db{isubj,1} '_' db{isubj,2} '_' svnm,'.mat'],'TRstruct', 'behavStruct','-v7.3')
    end
    clear at* behavStruct TRstruct

end
