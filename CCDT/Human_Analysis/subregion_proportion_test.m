% test proportions of selected features in each subregion
pretrial = 1; %=1 for pretrial features, =0 for intratrial features
featType = 1; % =0 for qexp & pow, =1 for qexp, =2 for pow
if pretrial == 1
%     load CCDTnodeAnalyze_PLV_80_pBP_qCAR_pre35_102323
%     load CCDTnodeAnalyze_PLV_80_pCAR_qCAR_pre40_hg70110_031324
%     load CCDTnodeAnalyze_PLV_80_pBP_qCAR_pre40_hg70110_031024 
%     load CCDTnodeAnalyze_PLV_100_pCAR_qCAR_pre_hg70110 %100 dat processed by me
%     load CCDTnodeAnalyze_PLV_100_pBP_qCAR_sIall_500mspre_origdat %orig data from vivek
%     load CCDTnodeAnalyze_PLV_100_pBP_qCAR_090324
%     load CCDTnodeAnalyze_PLV_100_pCAR_qCAR_090924_hg70110
%     load CCDTnodeAnalyze_PLV_80_pCAR_qCAR_hg70110_30_090924 
%     load CCDTnodeAnalyze_PLV_100_pCAR_qCAR_hg70110_090924.mat % 100dat
    load CCDTnodeAnalyze_PLV_80_pCAR_plvCAR_hg70110_30_031125 %plv instead of qexp
%     Qpolarity = 1; %1 = CAR, 2 = BP
%     Ppolarity = 1; %1 = CAR, 2 = BP
    polarity = 1; %both CAR
else
%     load CCDTnodeAnalyze_PLV_80_pBP_qCAR_intra30_102523
%     load CCDTnodeAnalyze_PLV_80_pCAR_qCAR_hg70-110_intra40_050624
%     load CCDTnodeAnalyze_PLV_80_pBP_qCAR_hg70-110_intra35_052224
%     load CCDTnodeAnalyze_PLV_80_pCAR_qCAR_hg70110_intra_30_090924
    load CCDTnodeAnalyze_PLV_100_pCAR_qCAR_hg70110_intra_090924.mat % 100dat
    Qpolarity = 1; %1 = CAR, 2 = BP
    Ppolarity = 1; %1 = CAR, 2 = BP
    polarity = 1;
end
load patient_loc_120623
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
            subreg_name_mat = repmat({EVE_roi_WM_csv{i}(1:end-4)},length(EVE_roi_WM_csv),1); %have to do end-8 to avoid inconsistent # underscores
        end
        X = find(cellfun(@startsWith,EVE_roi_WM_csv,subreg_name_mat));
        Ridx = [Ridx; X(find(X~=i))];
        if length(find(X~=i)) >1
            t = 1;
        end
    end
end

db = CCDTdatabase;
Nsubj = height(db);
Nregions = 8;
%
all_gm_wm_rois = []; all_wm_subrois = []; all_gm_subrois = [];
sig_chs_all_gm_wm_roi = []; sig_chs_all_wm_subrois = []; sig_chs_all_gm_subrois = [];
all_sig_chs_all_gm_wm_roi = []; all_sig_chs_all_wm_subrois = []; all_sig_chs_all_gm_subrois = [];

for isubj = 1:Nsubj
    if featType == 0
        sig_chs = [gchQ(ssIDq==isubj); gchP(ssIDp==isubj)];
    elseif featType == 1
        sig_chs = [gchQ(ssIDq==isubj)];
    elseif featType == 2
        sig_chs = [gchP(ssIDp==isubj)];
    end
    sig_chs_idx = zeros(size(sig_chs)); %index of significant channels
    %patient_loc(polarity).session(isubj).type~=0 fixes issue w/ subj 18
    ch_names = patient_loc(polarity).session(isubj).names(patient_loc(polarity).session(isubj).type~=0,1);
    for ch = 1:length(sig_chs)
        chidx = find(strcmp(ch_names,sig_chs(ch)));
        if ~isempty(chidx)
            sig_chs_idx(ch) = chidx;
        end
    end
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
sig_chs_all_wm_subrois(sig_chs_all_wm_subrois==0) = length(EVE_roi_WM_csv);
sig_chs_all_gm_subrois(sig_chs_all_gm_subrois==0) = length(seg_ID_list);
all_sig_chs_all_wm_subrois(all_sig_chs_all_wm_subrois==0) = length(EVE_roi_WM_csv);%w/ repeats
all_sig_chs_all_gm_subrois(all_sig_chs_all_gm_subrois==0) = length(seg_ID_list);%w/ repeats

%add ALIC/PLIC to TMwm region
% ALIC_PLIC_idx = all_wm_subrois==33 | all_wm_subrois==34 | all_wm_subrois==93 | all_wm_subrois==94;
% all_gm_wm_rois(ALIC_PLIC_idx) = 4;
% ALIC_PLIC_idx = sig_chs_all_wm_subrois==33 | sig_chs_all_wm_subrois==34 | sig_chs_all_wm_subrois==93 | sig_chs_all_wm_subrois==94;
% sig_chs_all_gm_wm_roi(ALIC_PLIC_idx) = 4;

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

% compare eve labels to channel count and aggregate roi assignment
% X = EVE_roi_WM_csv;
% for i = 1:length(EVE_roi_WM_csv)
%     X{i,2} = sum(all_wm_subrois==i);
%     X{i,3} = mean(all_gm_wm_rois(all_wm_subrois==i));
% end
% %%
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

%%
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

%% add percentages
for reg = 1:size(unique_subregions,2)
    format short
    % percent of all contacts in region
    unique_subregions{5,reg} = round(100*(unique_subregions{2,reg}/sum(unique_subregions{2,reg})),1);
    % percent of selected contacts in region
    unique_subregions{6,reg} = round(100*(unique_subregions{3,reg}/sum(unique_subregions{3,reg})),1);
    % percent of subregion selected 
    unique_subregions{7,reg} = round(100*(unique_subregions{3,reg}./unique_subregions{2,reg}),1);
end

%% TCWM pre vs post subregions
iReg = 4; %TCWM
for iSubreg = 1:6 %TCWM subregions
    X = [unique_subregions_pre{3,iReg}(iSubreg),unique_subregions_intra{3,iReg}(iSubreg)];
    N = [sum(unique_subregions_pre{3,iReg}),sum(unique_subregions_intra{3,iReg})];
    [h,p,chi2stat,df] = prop_test(X,N,false);
    disp(unique_subregions_pre{1,4}(iSubreg))
    p
end
% %pie chart
% color_assignment = [1,0,0; .3,.73,.05; .05,.4,.73; .9,.53,.1; 1,0,1; .27,.8,.85; .8,0,.5; 1,1,0]; %white background
% figure;
% set(gcf,'Color','white')
% p = pie(unique_subregions{3,4}');
% legend(["ACR","Postcentral","Precentral","PCR","Post. Thal. Rad.","SCR"])
% for i = 1:6
%     p(i*2-1).FaceColor = color_assignment(i,:);
% end
% title('Pretrial TCWM Subregion Selected Nodes')

% %line graph
% figure;
% set(gcf,'Color','white')
% title('TCWM Subregion Selected Nodes')
% hold on
% plot(unique_subregions_pre{3,4}([1,6,4,3,2,5]),'-o');% original order: "ACR","Postcentral","Precentral","PCR","Post. Thal. Rad.","SCR"
% plot(unique_subregions_intra{3,4}([1,6,4,3,2,5]),'-o');
% xticks(1:6)
% xticklabels(["ACR","SCR","PCR","Precentral","Postcentral","Post. Thal. Rad."])
% ylabel('% of TCWM Nodes from Subregion')
% legend('preparatory','anticipatory')
% xlim([0.5,6.5])

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
% ylim([0,50])
scatter([1,4],[47,47],'*','k')
% scatter([1,4],[79,79],'*','k') % stacked

% bar plot overlaid
% figure;
% set(gcf,'Color','white')
% b1 = bar(unique_subregions_pre{3,4}([1,6,4,3,2,5]),'r');
% b1.FaceAlpha = 0.5;
% hold on
% b2 = bar(unique_subregions_intra{3,4}([1,6,4,3,2,5]),.6,'b');
% b2.FaceAlpha = 0.5;
% xticks(1:6)
% xticklabels(["ACR","SCR","PCR","Precentral","Postcentral","Post. Thal. Rad."])
% ylabel('% of TCWM Nodes from Subregion')
% legend('preparatory','anticipatory')
% scatter([1,4],[47,47],'*','k')
% title('TCWM Subregion Selected Nodes')

%% qexp vs power selected features test
% % unique_subregions_q = unique_subregions;
% % unique_subregions_p = unique_subregions;
% qvp_test = [];
% for reg = 1:size(unique_subregions_q,2)
%     tot_region = sum(unique_subregions_q{2,reg});
%     tot_selected_q = sum(unique_subregions_q{3,reg}); 
%     tot_selected_p = sum(unique_subregions_p{3,reg}); 
%     X = [tot_selected_q, tot_selected_p];
%     N = [tot_region, tot_region];
%     [h,p,chi2stat,df] = prop_test(X,N,false);
%     qvp_test(reg,:) = [tot_selected_q/tot_region, tot_selected_p/tot_region, p, h];
% end

%% nonunique subregion plots
region_labels = {'Frontal GM', 'Temporoparietal GM', 'Paralimbic GM', ...
    'Thalamocortical WM','Frontal Association WM','Temporal Association WM',...
    'Paralimbic WM','Commissural Fibers'};
figure(1);
for i = 1:8 % # regions
    subplot(2,4,i)
    bar([nonunique_subregions{3,i}', unique_subregions{2,i}'])
    xticklabels(unique_subregions{1,i})
    title(region_labels{i})
    legend('selected features','non-selected features')
end


%% testing
% all_gm_wm_rois = [];
% all_all_gm_wm_rois = [];
% for isubj = 1:Nsubj
%     c_gm_wm_rois = [patient_loc(1).session(isubj).gm_wm_rois(patient_loc(1).session(isubj).type~=0);...
%         patient_loc(2).session(isubj).gm_wm_rois(patient_loc(2).session(isubj).type~=0)];
%     all_all_gm_wm_rois = [all_all_gm_wm_rois; c_gm_wm_rois];
%     c_all_chs = [patient_loc(1).session(isubj).names(patient_loc(1).session(isubj).type~=0);...
%         patient_loc(2).session(isubj).names(patient_loc(2).session(isubj).type~=0)];
%     [c_unique_chs, idx_unique_chs] = unique(c_all_chs);
%     c_gm_wm_rois = c_gm_wm_rois(idx_unique_chs);
%     disp(['# CAR chs: ' num2str(length(patient_loc(1).session(isubj).names(patient_loc(1).session(isubj).type~=0)))...
%         '  # BP chs: ' num2str(length(patient_loc(2).session(isubj).names(patient_loc(2).session(isubj).type~=0)))...
%         '  # unique chs: ' num2str(length(c_unique_chs))])
%     all_gm_wm_rois = [all_gm_wm_rois; c_gm_wm_rois];
% end
%% testing

% %% get subregions for all contacts
% Nregions = 8; 
% % frontal cortex = 1; temporoparietal cortex = 2; paralimbic GM = 3;
% % thalamocortical = 4; frontal assoc = 5; temporal assoc = 6; paralimbic = 7; commissural = 8;
% % load patient_loc_120623
% Nsubj = length(patient_loc(1).session);
% % gm seg ROI
% frontal_cortex = {'SFG', 'MFG', 'IFG', 'precentral', 'GRe', 'MOrG', 'POrG'};
% temporoparietal_cortex = {'MTG', 'ITG', 'STG', 'FuG', 'AnG', 'SMG', 'SPL','postcentral'};
% paralimbic_GM = {'Hippocampus','Amygdala','PHG','cingulate','LiG','Ent','insula'};
% seg_ID_list = [frontal_cortex,temporoparietal_cortex,paralimbic_GM,'n/a'];
% % get all ROI from EVE
% % EVE_roi_csv = readtable('EVE_all_roi.csv');
% EVE_roi_WM_csv = readtable('EVE_regions_WM.csv');
% % EVE_roi_list = EVE_roi_csv{2:end,2}; % extract all ROI
% % EVE_roi_list{end+1} = 'n/a';
% EVE_roi_WM_csv = EVE_roi_WM_csv{:,2};
% EVE_roi_WM_csv{end+1} = 'n/a';
% %concatenate all labels and rois across subjects
% all_wm_rois = [];
% all_gm_rois = [];
% all_gm_wm_rois = [];
% for isubj = 1:Nsubj
%     all_wm_rois = [all_wm_rois; int32(patient_loc(1).session(isubj).eve)];
%     all_gm_rois = [all_gm_rois; patient_loc(1).session(isubj).seg_ID(patient_loc(1).session(isubj).type~=0)];
%     all_gm_wm_rois = [all_gm_wm_rois; patient_loc(1).session(isubj).gm_wm_rois(patient_loc(1).session(isubj).type~=0)];
% end
% all_wm_rois(all_wm_rois==0) = length(EVE_roi_WM_csv);
% all_wm_labels = EVE_roi_WM_csv(all_wm_rois);
% all_gm_rois(all_gm_rois==0) = length(seg_ID_list);
% all_gm_labels = seg_ID_list(all_gm_rois);
% %%
% unique_subregions = cell(2,Nregions);
% for i = 1:3
%     subregions{i} = all_gm_labels(all_gm_wm_rois == i);
%     unique_subregions{1,i} = unique(all_gm_labels(all_gm_wm_rois == i));
%     Nsubregions = length(unique_subregions{1,i});
%     for j = 1:Nsubregions
%         unique_subregions{2,i}(j) = sum(ismember(subregions{i}, unique_subregions{1,i}(j)));
%     end
% end
% for i = 4:8
%     subregions{i} = all_wm_labels(all_gm_wm_rois == i);
%     unique_subregions{1,i} = unique(all_wm_labels(all_gm_wm_rois == i));
%     Nsubregions = length(unique_subregions{1,i});
%     for j = 1:Nsubregions
%         unique_subregions{2,i}(j) = sum(ismember(subregions{i}, unique_subregions{1,i}(j)));
%     end
% end
