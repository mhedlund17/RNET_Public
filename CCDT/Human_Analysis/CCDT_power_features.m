function CCDT_power_features
%   Compute narrowband and broadband power for each trial in a specified
%   perievent window. Generally follows analysis of Manning et al 2009.
%
%   DR 10/2019 & VB 5/2020 & MH 9/2024

% parameters
odir = ''; % output directory
svnm = ''; %name of output file
iev = 1; % event (1 = trial start cue, 2 = go cue, 3 = response)
win = [-500 0]; % perievent window (ms)
shiftDat = 0; %1 = create randomly shifted null data, 0= use regular data
fbands = [3 12; 12 30; 30 55; 70 110]; % frequency band ranges (Hz)
p.ddir = '/CCDT/eeg/'; % raw data directory
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
Nsubj = size(db,1); Nbands = size(fbands,1)+1;
load('patient_loc_101421.mat') %anatomical locations

% initiate structures to save
LTpow = cell(Nsubj,1); % learning statistic [m +/- 95% ci, p]
NFstruct = struct;

for isubj = stsubj:Nsubj
    disp([num2str(isubj) '/' num2str(Nsubj)]);
    p.subj = db{isubj,1}; q.subj = p.subj;
    p.sess = db{isubj,2}; q.sess = p.sess;
    q.ddir = [p.ddir(1:end-4) 'events/'];
    
    % load data
    [dat,Nsamp,fs,chn, gchlbl] = loadCCDTdata(p);

    % check contact labels to make sure they align w/ patient_loc file
    if length(gchlbl) == length(patient_loc(p.rrf).session(isubj).names)
        disp("gchs # = # chs in patient loc file")
    else
        disp("# gchs NOT the same as patient loc file")
        for ichlb = 1:length(patient_loc(p.rrf).session(isubj).names)
            for ichar = 1: length(patient_loc(p.rrf).session(isubj).names{ichlb})
                if gchlbl{ichlb}(ichar)~=patient_loc(p.rrf).session(isubj).names{ichlb}(ichar)
                    disp(['char mismatch in ch ' num2str(ichlb) '. channel removed.'])
                    gchlbl(ichlb) = [];
                    chn(ichlb) = [];
                    dat(:,ichlb) = [];
                end
            end
        end
        if length(gchlbl) == length(patient_loc(p.rrf).session(isubj).names)
            disp("gchs # = # chs in patient loc file")
        end
    end

    %remove contacts not labeled as gm or wm
    gchlbl(patient_loc(2).session(isubj).type==0',:)=[];
    chn(patient_loc(2).session(isubj).type==0)=[];
    dat(:,patient_loc(2).session(isubj).type==0) = [];
    Nch = size(dat,2);
    
    if shiftDat %randomly shift data in time series as null comparison
        ndat = zeros(size(dat,1),size(dat,2));
        for i=1:size(dat,2)
            ndat(:,i) = circshift(dat(:,i),randi(size(dat,1)));
        end
        dat = ndat;
    end
    
    
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
    RT = ccdt(:,5); % reaction time (ms)
    Ntrl = size(ccdt,1); % # trials
    
    % window indices for each trial relative to event/cue index
    swin = round(win(1)/1000*fs):round(win(2)/1000*fs);
    Nsamp = length(swin);
    indmat = ccdt(:,iev)*ones(1,Nsamp) + ones(Ntrl,1)*swin;
    
    % spectral power
    cPOW = zeros(Nch,Ntrl,Nbands);
    fb = cwtfilterbank('SignalLength', length(dat(:,1)), 'Wavelet', 'amor', 'SamplingFrequency',fs,'FrequencyLimits',[2 150], 'VoicesperOctave', 4);
    freq = centerFrequencies(fb);
    Nfreq = length(freq);
    for ii = 1:Nch
        datGPU = (dat(:,ii));
        P = cwt(datGPU,'FilterBank',fb);
        P = abs(P); %  power spectra
        P = P.^(1/4); % transform chi-square distributed values to normality (Hawkins and Wixley 1986)
        P = (P-mean(P(:)))./std(P(:)); % normalize to account for different electrode impedances

        % window spectra by trial
        Pwin = zeros(Nfreq,Ntrl);
        for jj = 1:Nfreq
            cP = P(jj,:);
            cPwin = cP(indmat);
            Pwin(jj,:) = mean(cPwin,2); % mean power in window
        end
        
        % narrowband power - for 1st 3 fbands
        for jj = 1:Nbands-1
            cbp = mean(Pwin(freq>=fbands(jj,1) & freq<=fbands(jj,2),:)); % mean power in band
            cPOW(ii,:,jj) = cbp;
        end
        
        % broadband power - for 4th (high gamma) fband
        for jj = 1:Ntrl
            b = robustfit(log10(freq),Pwin(:,jj));
            cy = b(1)+b(2)*log10(freq);
            cbbp = mean(cy);
            cPOW(ii,jj,Nbands) = cbbp;
        end
    end
    
    % flag trials with bad RTs
    inr = (RT<0 | RT>999); % incorrect responses, logical index
    disp([num2str(sum(inr)) ' trials omitted from analysis']);
    cTr = (~inr); %remove incorrect trials
    
    %regression - find features with a significant correlation with RT
    % includes all trials in this regression
    for ix=1:size(cPOW,1)
        for ij = 1:length(fbands)
            [bpow(ix,ij,:),statspow(ix,ij,:)] = robustfit(cPOW(ix,cTr,ij),RT(cTr));
            LTpow{isubj}(ix,ij,:) = [bpow(ix,ij,2) bpow(ix,ij,2)-1.96*statspow(ix,ij).se(2) bpow(ix,ij,2)+1.96*statspow(ix,ij).se(2) statspow(ix,ij).p(2)];
        end
    end
    
    for ij=1:length(fbands)
        NFstruct(isubj).fbands(ij).fbands = fbands(ij,:);
        NFstruct(isubj).fbands(ij).NFpow = cPOW(:,:,ij);
        NFstruct(isubj).gch = chn;
        NFstruct(isubj).gchlbl = gchlbl;
    end
    
    clear NFpow bpow stats*
end

if saveon
    save([odir svnm,'.mat'],'LT*','NFstruct')
end

end