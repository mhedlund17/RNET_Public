function CCDTpowRT2
%   Compute narrowband and broadband power for each trial in a specified
%   perievent window. Generally follows analysis of Manning et al 2009.
%
%   DR 10/2019 & VB 5/2020 & MH 9/2024

% parameters
p.ddir = '/mnt/sdb1/CCDT/eeg/'; % data directory
p.subj = []; % subject (leave empty to batch process all subjects in database)
p.sess = []; % session date (leave empty to chose automatically)
p.stime = []; % session time (leave empty for all session times)
p.rln = 1; % remove line noise? (0 or 1)
p.rrf = 1; % re-reference? (0 = none, 1 = common average re-reference, 2 = bipolar re-reference, 3 = laplacian re-reference)
p.outl = 1; % remove outlier (disconnected) channels? (0 or 1)
p.outlMetric = 'powGamma'; % outlier metric ('rms' = root-median-squared, 'powGamma' = high-gamma power)
p.outlThresh = 5; % outlier threshold (standard deviations)
sI = 1; %0 first task session only, 1 use all task sessions in same day
wavl = 'morl'; % wavelet
freq = logspace(log10(2),log10(150),50)'; % spectrum frequencies (Hz)
fbands = [3 12; 12 30; 30 55; 70 110]; % frequency range (Hz)
iev = 1 % event (1 = trial start cue, 2 = go cue, 3 = response)
win = [-500 0] % perievent window (ms)
shiftDat = 0; %create shifted null data
useMid = 0; %use only middle third of RT trials with linear regression
svnm = ''; %name of output file 
odir = ''; %output directory
stsubj = 1; % start subject
sigPThresh = 0.05; %significance threshold
saveon = 1; %save featuers

% load database
db = CCDTdatabase;
load('patient_loc_101421.mat') %anatomical locations
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

% loop through subjects
warning off
Nsubj = size(db,1); Nbands = size(fbands,1)+1; Nfreq = length(freq);
LTpow = cell(Nsubj,1); % learning statistic [m +/- 95% ci, p]
sP = cell(Nsubj,length(fbands),1); %index value within each subject good channels vector for most significant to least significant nodal Qexp RT prediction
sPval = cell(Nsubj,length(fbands),1);
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
    
    if shiftDat %randomly shift data
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
        if sI
            for ii = 2:length(ccdt)
                cccdt = ccdt{ii};
                cccdt(:,1:3) = cccdt(:,1:3) + sum(Nsamp(1:ii-1))*ones(size(cccdt,1),3); % event index accounting for concatenated data across session times
                nccdt = [nccdt; cccdt];
            end
        end
        ccdt = nccdt;
    end
    cDT = ccdt(:,4); % delay time (ms)
    cRT = ccdt(:,5); % reaction time (ms)
    Ntrl = size(ccdt,1);
    
    % window indices
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

        % window spectra
        Pwin = zeros(Nfreq,Ntrl);
        for jj = 1:Nfreq
            cP = P(jj,:);
            cPwin = cP(indmat);
            Pwin(jj,:) = mean(cPwin,2); % mean power in window
        end
        
        % narrowband power
        for jj = 1:Nbands-1
            cbp = mean(Pwin(freq>=fbands(jj,1) & freq<=fbands(jj,2),:)); % mean power in band
            cPOW(ii,:,jj) = cbp;
        end
        
        % broadband power
        for jj = 1:Ntrl
            b = robustfit(log10(freq),Pwin(:,jj));
            cy = b(1)+b(2)*log10(freq);
            cbbp = mean(cy);
            cPOW(ii,jj,Nbands) = cbbp;
        end
    end
    
    % flag trials with bad RTs
    inr = (cRT<0 | cRT>999); % incorrect responses, logical index
    disp([num2str(sum(inr)) ' trials omitted from analysis']);
    cTr = (~inr); %remove incorrect trials
    
    %regression
    for ix=1:size(cPOW,1)
        for ij = 1:length(fbands)
            [bpow(ix,ij,:),statspow(ix,ij,:)] = robustfit(cPOW(ix,cTr,ij),cRT(cTr));
            LTpow{isubj}(ix,ij,:) = [bpow(ix,ij,2) bpow(ix,ij,2)-1.96*statspow(ix,ij).se(2) bpow(ix,ij,2)+1.96*statspow(ix,ij).se(2) statspow(ix,ij).p(2)];
        end
    end
    for ij=1:length(fbands)
        [sPval{isubj}(ij,:),sP{isubj}(ij,:)] = sort(LTpow{isubj}(:,ij,4), 'ascend');
    end
    
    for ij=1:length(fbands)
        try
            infNodes = chn(LTpow{isubj}(:,ij,4)<=sigPThresh);
        catch
            infNodes = chn(sP{isubj}(ij,1));
            disp(['Freq ' num2str(ij) ' No sig nodes'])
        end
        NFstruct(isubj).fbands(ij).fbands = fbands(ij,:);
        NFstruct(isubj).fbands(ij).NFpow = cPOW(:,:,ij);
        NFstruct(isubj).gch = chn;
        NFstruct(isubj).gchlbl = gchlbl;
        NFstruct(isubj).fbands(ij).nodes = infNodes;
    end
    
    clear NFpow bpow stats* cN
end

if saveon
    save([odir svnm,'.mat'],'LT*','sP*','NFstruct')
end

warning on
end
