function [dat,Nsamp,fs,chn,chnm] = loadCCDTdata(varargin)
%   Load data recorded on Natus acquisition system in the EMU during CCDT
%   task. 
%
%   DR 03/2020

% parameters
ddir = []; % data directory
subj = []; % subject
sess = []; % session date
stime = []; % session time (leave empty for all session times)
outlThresh = 5; % outlier threshold (standard deviations)
if nargin
    v2struct(varargin{1});
end

% load data
p.ddir = ddir; p.subj = subj; p.sess = sess; p.stime = stime; p.verbose = 0;
p.reref = 'noreref'; % always use 'noreref' for CCDT analysis
[dat,Nsamp,fs,chn,chnm] = loadData_natus(p);

% remove non-sEEG (DC, EKG, EOG, EEG) channels
inseeg = find(strncmp(chnm,'DC',2)|strncmp(chnm,'EKG',3)|strcmp(chnm,'LOC')|strcmp(chnm,'ROC')|...
              strncmp(chnm,'C',1)|strncmp(chnm,'F',1)|strncmp(chnm,'P',1)|...
              strncmp(chnm,'O',1)|strncmp(chnm,'T',1)|strncmp(chnm,'FP',2));
if ~isempty(inseeg)
    dat(:,inseeg) = [];
    chn(inseeg) = [];
    disp([num2str(length(inseeg)) ' non-sEEG channels removed']);
    disp(chnm(inseeg)');
    chnm(inseeg) = [];
end

% remove null channels (no signal)
dat = dat - ones(size(dat,1),1)*median(dat,1); % remove dc offset
sigmaN = median(abs(dat)/.6745); % standard deviation of background noise
inull = find(sigmaN==0);
if ~isempty(inull)
    dat(:,inull) = [];
    chn(inull) = [];
    chnm(inull) = [];
    disp([num2str(length(inull)) ' null channels removed']);
end

Nch = size(dat,2); % # channels

% remove line noise
for ii = 1:Nch
    ln = mtmlinenoise(dat(:,ii),3,fs,fs,[60 120 180]);
    dat(:,ii) = dat(:,ii) - ln;
end

% remove outlier channels
s = bandpower(dat,fs,[80 200]); % high-gamma power (HFA activity)
zs = 0.6745*(s-median(s))/mad(s,1); % modified zscore of outlier statistic (Iglewicz and Hoaglin 1993)
ioutlier = find(abs(zs) > outlThresh); % Iglewicz and Hoaglin recommended 3.5 standard deviations - we chose 5 SD
if ~isempty(ioutlier)
    dat(:,ioutlier) = [];
    chn(ioutlier) = [];
    chnm(ioutlier) = [];
    disp([num2str(length(ioutlier)) ' outlier channels removed']);
end

% common average re-reference
dat = dat - mean(dat,2)*ones(1,size(dat,2));
