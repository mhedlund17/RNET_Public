% proportion test for selected features in each subregion

% set up
featType = 1; % =0 for qexp & pow, =1 for qexp, =2 for pow
load("preparatory file name") %load file saved in 2nd section of CCDTanalyze (with ssID and gch variables)
polarity = 1; %used to specify anatomical labels based on common average rereferencing in patient_loc file.
load patient_loc_120623 %anatomical localizations
db = CCDTdatabase; %subject database
Nsubj = height(db);
Nregions = 8;

% gm seg ROI
frontal_cortex = {'SFG', 'MFG', 'IFG', 'precentral', 'GRe', 'MOrG', 'POrG'};
temporoparietal_cortex = {'MTG', 'ITG', 'STG', 'FuG', 'AnG', 'SMG', 'SPL','postcentral'};
paralimbic_GM = {'Hippocampus','Amygdala','PHG','cingulate','LiG','Ent','insula'};
seg_ID_list = [frontal_cortex,temporoparietal_cortex,paralimbic_GM,'n/a'];

% wm ROI from EVE
EVE_roi_WM_csv = readtable('EVE_regions_WM.csv');
EVE_roi_WM_csv = EVE_roi_WM_csv{:,2};
EVE_roi_WM_csv{end+1} = 'n/a';

%get left and right subroi indices
Lidx = []; Ridx = [];
for i = 1:length(EVE_roi_WM_csv)-1
    if strcmp(EVE_roi_WM_csv{i}(end-3:end),'left')
        Lidx = [Lidx; i];
        if i ~= 62 %PONS creates error
            subreg_name_mat = repmat({EVE_roi_WM_csv{i}(1:end-8)},length(EVE_roi_WM_csv),1); %have to do end-8 to avoid inconsistent # underscores
        else
            subreg_name_mat = repmat({EVE_roi_WM_csv{i}(1:end-4)},length(EVE_roi_WM_csv),1); 
        end
        X = find(cellfun(@startsWith,EVE_roi_WM_csv,subreg_name_mat));
        Ridx = [Ridx; X(find(X~=i))];
        if length(find(X~=i)) >1
            t = 1;
        end
    end
end

%initiate empty variables
all_gm_wm_rois = []; all_wm_subrois = []; all_gm_subrois = [];
sig_chs_all_gm_wm_roi = []; sig_chs_all_wm_subrois = []; sig_chs_all_gm_subrois = [];
all_sig_chs_all_gm_wm_roi = []; all_sig_chs_all_wm_subrois = []; all_sig_chs_all_gm_subrois = [];

for isubj = 1:Nsubj
    if featType == 0 %power & qexp features
        sig_chs = [gchQ(ssIDq==isubj); gchP(ssIDp==isubj)];
    elseif featType == 1 %just qexp features
        sig_chs = [gchQ(ssIDq==isubj)];
    elseif featType == 2 %just power features
        sig_chs = [gchP(ssIDp==isubj)];
    end

    %match selected channel name to index in patient_loc file
    sig_chs_idx = zeros(size(sig_chs)); %index of significant channels
    %patient_loc(polarity).session(isubj).type~=0 fixes issue w/ subj 18 in patient_loc file
    ch_names = patient_loc(polarity).session(isubj).names(patient_loc(polarity).session(isubj).type~=0,1);
    for ch = 1:length(sig_chs)
        chidx = find(strcmp(ch_names,sig_chs(ch)));
        if ~isempty(chidx)
            sig_chs_idx(ch) = chidx;
        end
    end
    %channels names may be repeated if it was selected in multiple frequency bands
    sig_chs_idx_all = sig_chs_idx; % keep repeats to count all selected features per region
    sig_chs_idx = unique(sig_chs_idx); %remove repeats to count # nodes with selected features per region
    sig_chs_idx(sig_chs_idx==0) = []; %remove invalid channels
    sig_chs_idx_all(sig_chs_idx_all==0) = []; %remove invalid channels

    % get current subj rois
    c_gm_wm_rois = patient_loc(polarity).session(isubj).gm_wm_rois(patient_loc(polarity).session(isubj).type~=0);
    c_wm_subrois = int32(patient_loc(polarity).session(isubj).eve);
    c_gm_subrois = patient_loc(polarity).session(isubj).seg_ID(patient_loc(polarity).session(isubj).type~=0);
    
    %get composite regions
    all_gm_wm_rois = [all_gm_wm_rois; c_gm_wm_rois];
    sig_chs_all_gm_wm_roi = [sig_chs_all_gm_wm_roi; c_gm_wm_rois(sig_chs_idx)];
    all_sig_chs_all_gm_wm_roi = [all_sig_chs_all_gm_wm_roi; c_gm_wm_rois(sig_chs_idx_all)]; %w/ repeats

    %get significant ch's subregions for each composite region
    all_wm_subrois = [all_wm_subrois; c_wm_subrois];
    all_gm_subrois = [all_gm_subrois; c_gm_subrois];
    sig_chs_all_wm_subrois = [sig_chs_all_wm_subrois; c_wm_subrois(sig_chs_idx)];
    sig_chs_all_gm_subrois = [sig_chs_all_gm_subrois; c_gm_subrois(sig_chs_idx)];
    all_sig_chs_all_wm_subrois = [all_sig_chs_all_wm_subrois; c_wm_subrois(sig_chs_idx_all)]; %w/ repeats
    all_sig_chs_all_gm_subrois = [all_sig_chs_all_gm_subrois; c_gm_subrois(sig_chs_idx_all)]; %w/ repeats
end

%replace zeros with end indices for 'n/a' label
all_wm_subrois(all_wm_subrois==0) = length(EVE_roi_WM_csv);
all_gm_subrois(all_gm_subrois==0) = length(seg_ID_list);
sig_chs_all_wm_subrois(sig_chs_all_wm_subrois==0) = length(EVE_roi_WM_csv); %w/o repeats
sig_chs_all_gm_subrois(sig_chs_all_gm_subrois==0) = length(seg_ID_list); %w/o repeats
all_sig_chs_all_wm_subrois(all_sig_chs_all_wm_subrois==0) = length(EVE_roi_WM_csv);%w/ repeats
all_sig_chs_all_gm_subrois(all_sig_chs_all_gm_subrois==0) = length(seg_ID_list);%w/ repeats

% replace right eve ROI with left 
for i = 1:length(Ridx)
    all_wm_subrois(all_wm_subrois==Ridx(i)) = Lidx(i);
    sig_chs_all_wm_subrois(sig_chs_all_wm_subrois==Ridx(i)) = Lidx(i);
    all_sig_chs_all_wm_subrois(all_sig_chs_all_wm_subrois==Ridx(i)) = Lidx(i); %w/ repeats
end

%get labels for subrois
all_wm_subrois_labels = EVE_roi_WM_csv(all_wm_subrois);
all_gm_subrois_labels = seg_ID_list(all_gm_subrois);
sig_chs_all_wm_subrois_labels = EVE_roi_WM_csv(sig_chs_all_wm_subrois);
sig_chs_all_gm_subrois_labels = seg_ID_list(sig_chs_all_gm_subrois);
all_sig_chs_all_wm_subrois_labels = EVE_roi_WM_csv(all_sig_chs_all_wm_subrois); %w/ repeats
all_sig_chs_all_gm_subrois_labels = seg_ID_list(all_sig_chs_all_gm_subrois); %w/ repeats

% tally up contacts in each subregion
unique_subregions = cell(2,Nregions);
nonunique_subregions = cell(2,Nregions);
for i = 1:3 % gm rois
    all_subregions{i} = all_gm_subrois_labels(all_gm_wm_rois == i); %subregion label in each gm_wm_roi for all chs
    sig_chs_subregions{i} = sig_chs_all_gm_subrois_labels(sig_chs_all_gm_wm_roi == i); %subregion label in each gm_wm_roi for sig chs
    all_sig_chs_subregions{i} = all_sig_chs_all_gm_subrois_labels(all_sig_chs_all_gm_wm_roi == i); % w/ repeats
    unique_subregions{1,i} = unique(all_gm_subrois_labels(all_gm_wm_rois == i)); %name of subregions
    Nsubregions = length(unique_subregions{1,i});
    for j = 1:Nsubregions
        unique_subregions{2,i}(j) = sum(ismember(all_subregions{i}, unique_subregions{1,i}(j))); % total # chs w each subregion label
        unique_subregions{3,i}(j) = sum(ismember(sig_chs_subregions{i},unique_subregions{1,i}(j))); % # sig chs w each subregion label
        nonunique_subregions{3,i}(j) = sum(ismember(all_sig_chs_subregions{i},unique_subregions{1,i}(j))); % # sig chs w each subregion label
    end
end
for i = 4:8 % wm rois
    all_subregions{i} = all_wm_subrois_labels(all_gm_wm_rois == i); %subregion label in each gm_wm_roi for all chs
    sig_chs_subregions{i} = sig_chs_all_wm_subrois_labels(sig_chs_all_gm_wm_roi == i); %subregion label in each gm_wm_roi for sig chs
    all_sig_chs_subregions{i} = all_sig_chs_all_wm_subrois_labels(all_sig_chs_all_gm_wm_roi == i); % w/ repeats
    unique_subregions{1,i} = unique(all_wm_subrois_labels(all_gm_wm_rois == i));
    Nsubregions = length(unique_subregions{1,i});
    for j = 1:Nsubregions
        unique_subregions{2,i}(j) = sum(ismember(all_subregions{i}, unique_subregions{1,i}(j))); % total # chs w each subregion label
        unique_subregions{3,i}(j) = sum(ismember(sig_chs_subregions{i},unique_subregions{1,i}(j))); % # sig chs w each subregion label
        nonunique_subregions{3,i}(j) = sum(ismember(all_sig_chs_subregions{i},unique_subregions{1,i}(j))); % # sig chs w each subregion label
    end
end

% statistical testing
for reg = 1:size(unique_subregions,2)
    Nsubregions = length(unique_subregions{1,reg});
    tot_region = sum(unique_subregions{2,reg});
    tot_selected = sum(unique_subregions{3,reg});
    N = [tot_selected,tot_region];
    for subreg = 1:Nsubregions
        X = [unique_subregions{3,reg}(subreg),unique_subregions{2,reg}(subreg)];
        [h,p,chi2stat,df] = prop_test(X,N,false);
        unique_subregions{4,reg}(:,subreg) = [h;p];
    end
end

% add percentages
for reg = 1:size(unique_subregions,2)
    format short
    % percent of all contacts in region
    unique_subregions{5,reg} = round(100*(unique_subregions{2,reg}/sum(unique_subregions{2,reg})),1);
    % percent of selected contacts in region
    unique_subregions{6,reg} = round(100*(unique_subregions{3,reg}/sum(unique_subregions{3,reg})),1);
    % percent of subregion selected 
    unique_subregions{7,reg} = round(100*(unique_subregions{3,reg}./unique_subregions{2,reg}),1);
end

%% test if proportion of TCWM selected features are significantly different in prep vs antic periods
%run script above for preparatory and anticipatory features separately and 
% use lines below to keep both in the workspace for this comparison.

% unique_subregions_prep = unique_subregions;
% unique_subregions_antic = unique_subregions;

iReg = 4; %TCWM
for iSubreg = 1:6 %TCWM subregions
    X = [unique_subregions_prep{3,iReg}(iSubreg),unique_subregions_antic{3,iReg}(iSubreg)];
    N = [sum(unique_subregions_prep{3,iReg}),sum(unique_subregions_antic{3,iReg})];
    [h,p,chi2stat,df] = prop_test(X,N,false);
    disp(unique_subregions_prep{1,4}(iSubreg))
    p
end

%bar plot
figure;
set(gcf,'Color','white')
title('TCWM Subregion Selected Nodes')
hold on
bar([unique_subregions_pre{3,4}([1,6,4,3,2,5])',unique_subregions_intra{3,4}([1,6,4,3,2,5])']);% original order: "ACR","Postcentral","Precentral","PCR","Post. Thal. Rad.","SCR"
xticks(1:6)
xticklabels(["ACR","SCR","PCR","Precentral","Postcentral","Post. Thal. Rad."])
ylabel('% of TCWM Nodes from Subregion')
legend('preparatory','anticipatory')
xlim([0.5,6.5])
scatter([1,4],[47,47],'*','k')
