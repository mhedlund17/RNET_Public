%% CCDT_anatomical_analysis.m
% John Bernabei
% Dec 2020
% Center for Neuroengineering & Therapeutics
% 
% MH 08/23

%% 1. set up workspace
clear all

% feature data info
pdir = ''; % feature directory
pfname = ''; % power feature file name
gfname = ''; % graph communicability feature file name

% load feature selection info
sfeat_dir = ''; %feature selection directory
sfeat_fname = ''; %feature selection file name - files saved in CCDT_feature_select.m
sfeatInfo_dir = ''; %directory
sfeatInfo_fname = ''; %file saved in CCDT_feature_details.m
load([sfeatInfo_dir sfeatInfo_fname])
percThresh = 0.30; % percentage, as a decimal 
% percent threshold: features with a significant regression slope in
% more than this percentage of bootstrap iterations will be included in
% the selected feature space. The threshold used for analysis in paper was
% chosen by evaluating SVM performance in predicting fast vs slow trials
% for 20 thresholds between 0 and 100%. 

%parameters for heatmap figure output
addAsterisk = 1; %add asterisk to heatmap figure
save_heatmap_data = 0; %save heatmap figure info
savedir = ''; %output directory 
savefname = ''; % output file name

% get all ROI from EVE white matter atlas
EVE_roi_csv = readtable('EVE_all_roi.csv');
EVE_roi_WM_csv = readtable('EVE_regions_WM.csv');
EVE_roi_list = EVE_roi_csv{:,2}; % extract all ROI

% session list - 27 sessions representing 24 patients with SEEG
db = CCDTdatabase;
Nsubj = height(db);

% define colors for plots
gm_color_1 = [0.02, 0.62, 0.584]; %teal
gm_color_2 = [0.173, 0.839, 0.796];
gm_color_3 = [0.486, 0.98, 0.949];
wm_color_1 = [0.278, 0.059, 0.38];
wm_color_2 = [0.439, 0.024, 0.631]; %purple
wm_color_3 = [0.608, 0.259, 0.769];
wm_color_4 = [0.855, 0.545, 1];
wm_color_5 = [0.937, 0.804, 1];
wm_color_6 = [0.973, 0.91, 1];
barcolors = [gm_color_1;gm_color_2;gm_color_3;wm_color_2;wm_color_3;wm_color_4;wm_color_5;wm_color_6]; %gm/wm each diff shades of same color

%% 2. create EVE roi list bilateral
% we're taking all regions besides 'background' and combining R + L
EVE_roi_list_bilat{1,1} = 'Background';
for i = 2:89
    split_regions = split(EVE_roi_list{i},'_left'); % some string parsing
    EVE_roi_list_bilat(i,1) = split_regions(1);
end
    
% list which ROI of the Hammer smith atlas belong to which lobe (frontal vs temporal)
temporal_list = [1:16,30:31,82:83];
frontal_list = [20:21, 24:29, 50:59, 68:73, 76:81];

% list top segmentation regions of GM to use based on prevalence
seg_gm_list = {'MFG','MTG','ITG','fusiform','Hippocampus'};

% List EVE regions for each composite region = white matter
thalamo_cortical = [36,37,38,6,7,35,96,97,98,68,69,95]; 
frontal_association = [3,4,5,42,43,46,21,24,65,66,67,102,103,106,82,85];
temporal_association = [20,19,18,45,44,12,8,1,81,70,79,105,104,74,70,64];
paralimbic_wm = [41,39,40,11,47,101,100,99,73,107];
commissural = [51,52,53,111,112,113];

% list segmentation labels for each composite region = grey matter
frontal_cortex = {'SFG', 'MFG', 'IFG', 'precentral', 'GRe', 'MOrG', 'POrG'};
temporoparietal_cortex = {'MTG', 'ITG', 'STG', 'FuG', 'AnG', 'SMG', 'SPL','postcentral'};
paralimbic_GM = {'Hippocampus','Amygdala','PHG','cingulate','LiG','Ent','insula'};

% List EVE regions for each composite region = grey matter
frontal_cortex_eve = [3,4,5,6,22,24,91,92,93,94,109,112];
temporoparietal_cortex_eve = [1,7,8,12,18,19,20,23,89,95,96,100,106,107,108,111];
paralimbic_gm_eve = [2,11,13,17,25,26,27,90,99,101,105,113,114,115];

% assign the individual region lists into a structure
% this is the general ROI order we keep throughout
all_sub_roi(1).data = frontal_cortex_eve;
all_sub_roi(2).data = temporoparietal_cortex_eve;
all_sub_roi(3).data = paralimbic_gm_eve;
all_sub_roi(4).data = thalamo_cortical;
all_sub_roi(5).data = frontal_association;
all_sub_roi(6).data = temporal_association;
all_sub_roi(7).data = paralimbic_wm;
all_sub_roi(8).data = commissural;

%% 3. Electrode processing and localization for all patients
% *** not necessary to do if patient_loc file already loaded ***

polarity = 1; % Common Average Rereferencing done
% do separate localizations for monopolar & bipolar montage
for featureType = ["qexp", "pow"] 
    
    % need different labels and thus localization for power vs qexp
    if strcmp(featureType,"qexp")
        % getting channel labels from graphRT structure
        load([pdir gfname])        
        
        % set up empty vectors for containing the following:
        quad_nodes_q = []; % whether each monopolar node is in GM or WM
        quad_lobes_q = []; % which lobe each node is in (frontal vs temporal)
        quad_seg_q = []; % which segmentation ROI each node is in
        quad_wm_q = []; % which EVE WM structure each node is in
    else 
        % getting channel labels from powRT structure
        load([pdir pfname])
        
        % set up empty vectors for containing the following:
        quad_nodes_p = []; % whether each bipolar node is in GM or WM
        quad_lobes_p = []; % which lobe each node is in (frontal vs temporal)
        quad_seg_p = []; % which segmentation ROI each node is in
        quad_wm_p = []; % which EVE WM structure each node is in
    end 
    
    % T1 segmentation ROI structure containing all ROI across patients
    all_T1_roi = [];
 
    % do basic localization of all nodes
    for session = [1:27]
        disp([char(featureType) ' Subject: ' num2str(session)])
       
        % get patient
        patient = db{session,1};
        
        % clear out variables
        clear this_T1_regions type_ind elec_mni_coords seg_roi seg_roi1 seg_roi2 this_yeo_network yeo_network
        clear JHU_wm_list node_type lobe_type seg_type wm_type lobe_list gm_seg_inds this_wm_label
        
        % add these files to the github
        % Load electrode names & coordinates in MNI space from CSV file
        elec_label_mni_raw = readtable(sprintf('CCDTmni/new_%s_localization.csv',patient));
        
        % extract channel labels
        ch_names = NFstruct(session).gchlbl; 
        final_ch_names = ch_names;
        
        v = 1;
        % loop through channels and assign ROI & coords to parse ch
        % names to labels from localization
        for e = 1:length(ch_names)
            
            % extract electrode name
            elec1 = ch_names{e,1};
            
            % match up names
            coord_ind_1 = find(strcmp(elec_label_mni_raw{:,1},elec1));
            
            % extract MNI coordinates
            mni_1 = elec_label_mni_raw{coord_ind_1,3:5};
            
            % get the MNI coords and segmentation ROI
            try elec_mni_coords(v,:) = mni_1;
                seg_roi(v,1) = elec_label_mni_raw{coord_ind_1,2};
                yeo_network(v,1) = elec_label_mni_raw{coord_ind_1,7};
                v = v+1;
            catch ME
                % print if electrode is not available 
                % (segmentation error where a small number of contacts are unlocalized)
                fprintf('could not find electrode\n')
                final_ch_names(e,:) = [];
            end
        end
        
        % get number of electrodes
        num_elecs = size(elec_mni_coords,1);
        
        % determine electrode type: 0 = gm, 1 = wm, 2 = out of brain
        type_ind = zeros(num_elecs,1);
        
        % Correct ROI to go from blank to 'n/a'
        for e = 1:num_elecs
            if strcmp(seg_roi{e},'')
                this_T1_regions(e,1) = {'n/a'};
                this_yeo_network(e,1) = {'n/a'};
                type_ind(e) = 0;
            else
                this_T1_regions(e,1) = seg_roi(e);
                this_yeo_network(e,1) = yeo_network(e);
                type_ind(e) = 1; % this means it is in the brain
            end

            if strcmp(this_T1_regions{e,1},'Left Cerebral White Matter')||strcmp(this_T1_regions{e,1},'Right Cerebral White Matter')
                type_ind(e) = -1; % change to -1 if it is in white matter
            end
        end
        
        % assign GM segmentation indices (first three composite ROI)
        gm_seg_inds = zeros(length(this_T1_regions),1);
        for r_gm = 1:length(gm_seg_inds)
            this_region = this_T1_regions{r_gm}; % from seg localization
            
            if contains(this_T1_regions{r_gm},frontal_cortex)
                gm_seg_inds(r_gm) = 1; % frontal cortex = 1
                
            elseif contains(this_T1_regions{r_gm},temporoparietal_cortex)
                gm_seg_inds(r_gm) = 2; % temporoparietal cortex = 2
                
            elseif contains(this_T1_regions{r_gm},paralimbic_GM)
                gm_seg_inds(r_gm) = 3; % paralimbic GM = 3 
            end
        end
        
        % need individual label ID in case we do a sub-ROI analysis
        gm_seg_label_ID = zeros(length(this_T1_regions),1);
        for r_gm = 1:length(gm_seg_label_ID)
            this_region = this_T1_regions{r_gm};
            % frontal
            for r_atl = 1:length(frontal_cortex)
                if contains(this_T1_regions{r_gm},frontal_cortex{r_atl})
                    gm_seg_label_ID(r_gm) = r_atl;
                end
            end
            % temporoparietal
            for r_atl = 1:length(temporoparietal_cortex)
                if contains(this_T1_regions{r_gm},temporoparietal_cortex{r_atl})
                    gm_seg_label_ID(r_gm) = r_atl+length(frontal_cortex);
                end
            end
            % paralimbic
            for r_atl = 1:length(paralimbic_GM)
                if contains(this_T1_regions{r_gm},paralimbic_GM{r_atl})
                    gm_seg_label_ID(r_gm) = r_atl+length(frontal_cortex)+length(temporoparietal_cortex);
                end
            end 
        end
        
        % calculate number of GM & number of WM
        num_gm(session) = sum(type_ind==1); % total number of GM for this pt
        num_wm(session) = sum(type_ind==-1); % total number of WM for this pt
        num_ex_brain(session) = sum(type_ind==0); % total number of unlocalized
        
        % get which electrodes are usable (GM + WM)
        usable_inds = (type_ind==1)+(type_ind==-1);
        
        % JHU EVE atlas - GM+WM - localization 
        [mni_coords, node_roi, NN_flag] = nifti_values(elec_mni_coords(usable_inds==1,:),'JHU_MNI_SS_WMPM_Type-I_to_MNI_brain.nii');
        
        % JHU EVE atlas - WM only - localization
        [mni_wm_only, EVE_wm_only, NN_wm_only] = nifti_values(elec_mni_coords(usable_inds==1,:),'JHU_MNI_SS_WMPM_Type-III_to_MNI_brain.nii');

        % Hammer smith atlas
        [mni_lobar, lobar_roi, NN_lobar] = nifti_values(elec_mni_coords(usable_inds==1,:),'Hammers_mith_atlas_n30r83_SPM5.nii');
        
        % extract which lobe each lobe belongs to
        for l = 1:length(lobar_roi)
            if ismember(lobar_roi(l),temporal_list)
                lobe_list(l) = 2; % 2 = temporal lobe
            elseif ismember(lobar_roi(l),frontal_list)
                lobe_list(l) = 1; % 1 = frontal lobe
            else
               lobe_list(l) = 0; % 0 = other
            end
        end
        
        % extract which WM structure each EVE electrode belongs to
        for l = 1:length(EVE_wm_only)
            this_wm_roi = EVE_wm_only(l); % get the ROI
            try this_wm_label(l,1) = EVE_roi_WM_csv{EVE_wm_only(l),2};
            catch anyerror
                this_wm_label(l,1) = {'n/a'};
            end
            % check if it is in any of these sub - ROIs
            if ismember(this_wm_roi,thalamo_cortical)
                JHU_wm_list(l,1) = 1; 
            elseif ismember(this_wm_roi,frontal_association)
                JHU_wm_list(l,1) = 2;
            elseif ismember(this_wm_roi,temporal_association)
                JHU_wm_list(l,1) = 3;
            elseif ismember(this_wm_roi,paralimbic_wm) 
                JHU_wm_list(l,1) = 4;
            elseif ismember(this_wm_roi,commissural)
                JHU_wm_list(l,1) = 5;
            else 
                JHU_wm_list(l,1) = 0;
            end
        end

        % for all indices in brain, assign a WM ROI from EVE (GM regions should be 0)
        all_wm_roi = zeros(length(type_ind),1);
        all_wm_roi(find(usable_inds),1) = JHU_wm_list;
        
        % for all indices in brain, assign a GM/WM ROI from EVE
        all_node_roi = zeros(length(type_ind),1);
        all_node_roi(find(usable_inds),1) = node_roi;
        
        % for all indices in brain, assign a lobe
        all_lobe_list = zeros(length(type_ind),1);
        all_lobe_list(find(usable_inds),1) = lobe_list;
        
        % now we must make sure we only keep wm roi which are definitely wm
        % based on the segmentation
        all_non_wm = [find(type_ind==1);find(type_ind==0)];
        all_wm_roi(all_non_wm) = 0;
        
        % final gm/wm rois to use in paper
        ROI_labels = {'Frontal GM','Temporal GM','Paralimbic GM','Thalamocortical WM',...
            'Frontal Assoc WM','Temporal Assoc WM','Paralimbic WM','Commissural WM'};
        gm_wm_rois = gm_seg_inds; % first 3 rois same as gm_seg_inds
        for iwm = 1:5 %wm rois
            gm_wm_rois(all_wm_roi==iwm) = iwm+3;
        end
        gm_wm_labels = cell(length(gm_wm_rois),1);
        gm_wm_labels(gm_wm_rois~=0) = ROI_labels(gm_wm_rois(gm_wm_rois~=0));

        % create a patient localization structure
        patient_loc(polarity).session(session).coords = elec_mni_coords; % mni coordinates
        patient_loc(polarity).session(session).names = final_ch_names; % ch names
        patient_loc(polarity).session(session).type = type_ind; % wm vs gm vs outside brain
        patient_loc(polarity).session(session).roi = all_node_roi; % EVE GM & WM localization
        patient_loc(polarity).session(session).lobes = all_lobe_list; % lobe list from hammer-smith
        patient_loc(polarity).session(session).t1 = this_T1_regions; % t1 region from segmentation
        patient_loc(polarity).session(session).yeo = this_yeo_network; % yeo network from mni
        patient_loc(polarity).session(session).seg = gm_seg_inds; % GM segmentation inds ()1-3
        patient_loc(polarity).session(session).eve = EVE_wm_only'; % EVE WM structures
        patient_loc(polarity).session(session).seg_ID = gm_seg_label_ID; % GM sub-ROI
        patient_loc(polarity).session(session).wm = all_wm_roi; % WM structures
        patient_loc(polarity).session(session).wm_label = this_wm_label; % WM structures
        patient_loc(polarity).session(session).gm_wm_rois = gm_wm_rois; % GM + WM ROIs (used in final paper)
        patient_loc(polarity).session(session).gm_wm_labels = gm_wm_labels; % GM + WM ROI labels (used in final paper)
           
    end
end

%% 4. Do anatomical analysis
% scatter plot figure for comparing nonselected features to
% selected features in each anatomical region

seg_counts_nonsig = NaN*zeros(27,8);
patient_vec = [1:27];

for r = 1:8
    which_subregions = all_sub_roi(r).data;
    all_region_counts(r).data = zeros(length(which_subregions),2,8);
    seg_sig_coord_struct(r).data = [];
    seg_ns_coord_struct(r).data = [];
end

plot_x_label = {'Frontal','Tmp-parietal','Paralimbic',...
                'Thal.-cortical','Frontal Assoc.','Tmp.Assoc.',...
                'Paralimbic','Commissural'};

qexp_vs_pow = {'qexp','pow'};
threshold = [0.05, 1.1, 1.1];
color1 = [0, 0.4470, 0.7410];
color2 = [0.6350, 0.0780, 0.1840];
color5 = [103 55 155]/255;
color6 = [78 172 91]/255;
folder_list = {'p_005','all','min_max'};
freq_bands = {'Alpha-theta','Beta','Low-gamma','High-gamma'};

% get b and p values from feature selection
load([sfeat_dir,sfeat_fname])
for session = 1:height(sigPow)
    powB{session} = mean(sigPow{session,3},3);
    comB{session} = mean(sigCom{session,3},3);

    sigBpow = NaN*zeros(size(sigPow{session,3}));
    sigBpow(sigPow{session,1}) = sigPow{session,3}(sigPow{session,1});
    sigBcom = NaN*zeros(size(sigCom{session,3}));
    sigBcom(sigCom{session,1}) = sigCom{session,3}(sigCom{session,1});
    powBsig{session} = mean(sigBpow,3,'omitnan');
    comBsig{session} = mean(sigBcom,3,'omitnan');

    sigPpow = NaN*zeros(size(sigPow{session,2}));
    sigPpow(sigPow{session,1}) = sigPow{session,2}(sigPow{session,1});
    sigPcom = NaN*zeros(size(sigCom{session,2}));
    sigPcom(sigCom{session,1}) = sigCom{session,2}(sigCom{session,1});
    powP{session} = mean(sigPpow,3,'omitnan');
    comP{session} = mean(sigPcom,3,'omitnan');
end
         
z = 0;
for metric = 1:2 % do for qexp and pow
    clear all_pval EVE_all_freq
    
    %load feature data
    if metric==1
        load([pdir gfname])
    else 
        load([pdir pfname])
    end
     
    %initiate empty structures
    sig_ch_all = cell(Nsubj,1);
    sig_ch_names_all = cell(Nsubj,1);
    
    for freq = 1:4 %4 freq bands
        z = z+1;
        
        for r = 1:8
            seg_bval_ns(r).data = [];
        end
        
        % do for thresholds 0.05 and 0.1
         clear patient_roi
         clear mean_bval_lobe
         clear EVE_100_elec
             
         all_non_sig = [];
         pt_non_sig = [];
         seg_sig_coords = [];
         wm_sig_coords = [];
         seg_sig_all = [];
         wm_sig_all = [];
         background_sig_coords = [];
         background_pval_all = [];
         seg_pval_all = [];
         wm_pval_all = [];
         roi_sig_coords = [];
         wm_roi_sig_coords = [];
             
         for r = 1:9
            seg_EVE_roi(r).data = [];
            ns_seg_EVE_roi(r).data = [];
         end
             
         patient_roi = NaN*zeros(1,27);
         seg_bvals = NaN*zeros(27,3);
         seg_counts = NaN*zeros(27,9);
         EVE_100_elec = NaN*zeros(27,12);
         patient_roi_min = zeros(177,27);
         patient_roi_max = zeros(177,27);
            
        % do for sessions 1:27
        for session = patient_vec %session
            % Get all b and p values
            if strcmp(qexp_vs_pow{metric},'qexp')
                all_b_val = comB{session}(:,freq);
                all_p_val = comP{session}(:,freq);
                all_b_val_sig = comBsig{session}(:,freq);
                sigPerc = sum(sigCom{session,1}(:,freq,:),3)/size(sigCom{session,1},3);
            else
                all_b_val = powB{session}(:,freq); %avg b val avged over all iterations
                all_p_val = powP{session}(:,freq);
                all_b_val_sig = powBsig{session}(:,freq); %avg b val avged over iterations w p<.05
                sigPerc = sum(sigPow{session,1}(:,freq,:),3)/size(sigPow{session,1},3);
            end
                
            % get significant channels based on p value threshold
            [sig_ch,sigFreq] = find(sigPerc>percThresh);
            [non_sig_ch,nonsigFreq] = find(sigPerc<=percThresh);
            sig_ch_all{session} = [sig_ch_all{session}; sig_ch];
            ch_names = NFstruct(session).gchlbl;
                                
            % check if there are any significant channels
            if isempty(sig_ch)
                fprintf('skipping session b/c no significance\n')
                pt_non_sig = [pt_non_sig;mean(all_b_val)];
                all_non_sig = [all_non_sig;all_b_val];
            else
                % extract significant channel p values, p values, and names
                sig_b_val = all_b_val_sig(sig_ch);
                sig_p_val = all_p_val(sig_ch);
                sig_ch_names = ch_names(sig_ch);
                nonsig_ch_names = ch_names(non_sig_ch);
                sig_ch_names_all{session} = [sig_ch_names_all{session}; sig_ch_names];
                
                non_sig_b_val = all_b_val(non_sig_ch);
                pt_non_sig = [pt_non_sig;mean(non_sig_b_val,'omitnan')];
                all_non_sig = [all_non_sig;non_sig_b_val];
                
                patient_coord_mapping = zeros(length(sig_ch_names),1);
                ns_patient_coord_mapping = zeros(length(nonsig_ch_names),1);
                for ch = 1:length(sig_ch_names)
                    % get this individual channel name
                    this_ch = sig_ch_names{ch,:};
                    % get coordinate index
                    coord_index = find(strcmp(patient_loc(polarity).session(session).names(:,1),this_ch));
                    
                    try patient_coord_mapping(ch) = coord_index;
                    catch ME
                        %fprintf('cannot find electrode %s\n', this_ch);     
                    end
                end
                    
                for ch = 1:length(nonsig_ch_names)
                    % get this individual channel name
                    this_ch = nonsig_ch_names{ch,:};
                    % get coordinate index
                    coord_index = find(strcmp(patient_loc(polarity).session(session).names(:,1),this_ch));
                    
                    try ns_patient_coord_mapping(ch) = coord_index;
                    catch ME
                        %fprintf('cannot find electrode %s\n', this_ch);
                        
                    end
                end
                    
                [ind_max] =  (sig_b_val==max(sig_b_val));
                [ind_min] =  (sig_b_val==min(sig_b_val));
                    
                % need csv files for all, <0.05, and min/max
                    
                sig_roi = patient_loc(polarity).session(session).roi(patient_coord_mapping(patient_coord_mapping~=0));
                sig_type = patient_loc(polarity).session(session).type(patient_coord_mapping(patient_coord_mapping~=0));
                sig_lobe = patient_loc(polarity).session(session).lobes(patient_coord_mapping(patient_coord_mapping~=0));
                sig_seg = patient_loc(polarity).session(session).seg(patient_coord_mapping(patient_coord_mapping~=0));
                
                nonsig_seg = patient_loc(polarity).session(session).seg(ns_patient_coord_mapping(ns_patient_coord_mapping~=0));
                nonsig_wm = patient_loc(polarity).session(session).wm(ns_patient_coord_mapping(ns_patient_coord_mapping~=0));
                
                sig_wm = patient_loc(polarity).session(session).wm(patient_coord_mapping(patient_coord_mapping~=0));
                sig_mni = patient_loc(polarity).session(session).coords(patient_coord_mapping(patient_coord_mapping~=0),:);
                ns_mni = patient_loc(polarity).session(session).coords(ns_patient_coord_mapping(ns_patient_coord_mapping~=0),:);
                %sig_eve = patient_loc(metric).session(session).eve(patient_coord_mapping(patient_coord_mapping~=0),:);
                sig_seg_ID = patient_loc(polarity).session(session).seg_ID(patient_coord_mapping(patient_coord_mapping~=0),:);
                proc_pval = sig_p_val(patient_coord_mapping~=0);
                proc_bval = sig_b_val(patient_coord_mapping~=0);
                %ns_proc_pval = non_sig_p_val(patient_coord_mapping~=0);
                ns_proc_bval = non_sig_b_val;%(patient_coord_mapping~=0);
                    
                % do GM/WM analysis
                mean_bval(session,(2*(freq-1)+1)) = mean(proc_bval(sig_type==-1));
                mean_bval(session,(2*freq)) = mean(proc_bval(sig_type==1));
                    
                frac_sig(session,(2*(freq-1)+1)) = sum(sig_type==-1)./sum(patient_loc(polarity).session(session).type==-1);
                frac_sig(session,(2*freq)) = sum(sig_type==1)./sum(patient_loc(polarity).session(session).type==1);
                    
                % GM/WM + lobe analysis
                mean_bval_lobe(session,1) = mean(proc_bval([sig_type==-1]&[sig_lobe==1]));
                mean_bval_lobe(session,2) = mean(proc_bval([sig_type==1]&[sig_lobe==1]));
                mean_bval_lobe(session,3) = mean(proc_bval([sig_type==-1]&[sig_lobe==2]));
                mean_bval_lobe(session,4) = mean(proc_bval([sig_type==1]&[sig_lobe==2]));
                    
                d = 0;
                for r = [0 4 12 19 20 27 37 68 69 70 82 83]
                    d = d+1;
                    if r==0
                        background_vals = sig_b_val(sig_roi==r);
                        EVE_100_elec(session,d) = mean(background_vals,'omitnan');
                    else
                        left_side = sig_b_val(sig_roi==r);
                        right_side = sig_b_val(sig_roi==(r+88));
                        EVE_100_elec(session,d) = mean([left_side;right_side],'omitnan');
                    end
                end
                    
                % go through GM segmentation aggregate regions
                for r = 1:3
                    seg_bvals(session,r) = mean(sig_b_val(sig_seg==r),'omitnan');
                    seg_counts(session,r) = length(sig_b_val(sig_seg==r));
                    seg_counts_nonsig(session,r) = length(non_sig_b_val(nonsig_seg==r)); %added new
                    seg_bval_ns(r).data = [seg_bval_ns(r).data; non_sig_b_val(nonsig_seg==r)];
                    seg_EVE_roi(r).data = [seg_EVE_roi(r).data;sig_b_val(sig_seg==r)];
                    seg_sig_all = [seg_sig_all; sig_b_val(sig_seg==r)];
                    seg_pval_all = [seg_pval_all; sig_p_val(sig_seg==r)];
                    seg_sig_coords = [seg_sig_coords; sig_mni((sig_seg==r),:)];
                    roi_sig_coords = [roi_sig_coords; r*ones(size(sig_mni((sig_seg==r),:),1),1)];
                    
                    seg_sig_coord_struct(r).data =  [seg_sig_coord_struct(r).data; sig_mni((sig_seg==r),:)];
                    seg_ns_coord_struct(r).data =  [seg_ns_coord_struct(r).data; ns_mni((nonsig_seg==r),:)];
                    
                    ns_seg_EVE_roi(r).data = [ns_seg_EVE_roi(r).data;non_sig_b_val(nonsig_seg==r)];
                end
                    
                % go through WM EVE aggregate regions
                for r = 1:5
                    EVE_WM(session,r) = mean(sig_b_val(sig_wm==r),'omitnan');
                    seg_counts(session,(r+3)) = length(sig_b_val(sig_wm==r));
                    seg_counts_nonsig(session,(r+3)) = length(non_sig_b_val(nonsig_wm==r));
                    seg_bval_ns(r+3).data = [seg_bval_ns(r+3).data; non_sig_b_val(nonsig_wm==r)];
                    seg_EVE_roi(r+3).data = [seg_EVE_roi(r+3).data;sig_b_val(sig_wm==r)];
                    wm_sig_all = [wm_sig_all; sig_b_val(sig_wm==r)];
                    wm_pval_all = [wm_pval_all; sig_p_val(sig_wm==r)];
                    wm_sig_coords = [wm_sig_coords; sig_mni((sig_wm==r),:)];
                    wm_roi_sig_coords = [wm_roi_sig_coords; (r+3)*ones(size(sig_mni((sig_wm==r),:),1),1)];
                    
                    seg_sig_coord_struct(r+3).data = [seg_sig_coord_struct(r+3).data; sig_mni((sig_wm==r),:)];
                    seg_ns_coord_struct(r+3).data =  [seg_ns_coord_struct(r+3).data; ns_mni((nonsig_wm==r),:)];
                    
                    ns_seg_EVE_roi(r+3).data = [ns_seg_EVE_roi(r+3).data;non_sig_b_val(nonsig_wm==r)];
                end

                seg_counts(session,9) = length(sig_b_val(sig_type==0));
                
                pt_background(session,1) = mean(sig_b_val(sig_type==0),'omitnan');
                seg_EVE_roi(9).data = [seg_EVE_roi(9).data;sig_b_val(sig_type==0)];
                
                background_sig_coords = [background_sig_coords; sig_mni((sig_type==0),:)];
                background_pval_all = [background_pval_all;sig_p_val(sig_type==0)];
                
                % deterine how many nodes are in each roi
                for r = 0:176
                    row_ind = r+1;
                    patient_roi(row_ind,session) = sum(sig_roi==r);
                end
                    
            end
            
            all_vals_sig(z).seg = seg_sig_all;
            all_vals_sig(z).wm = wm_sig_all;
        end
        
        feat_nums = sum(seg_counts(:,1:9),2);
        feat_nums(isnan(feat_nums)) = 0;
        feat_gen(1:27,z) = feat_nums(1:Nsubj);
            
        all_non_sig_struct(z).data = all_non_sig;
        pt_non_sig_struct(z).data = pt_non_sig;
        
        EVE_all_freq(:,:,freq) = EVE_100_elec;
            
        for r = 1:8
            all_median_table(z,r) = nanmedian(seg_EVE_roi(r).data);  
            all_mean_table(z,r) = mean(seg_EVE_roi(r).data,'omitnan');
            all_median_table_nonsig(z,r) = nanmedian(all_non_sig_struct(z).data);
            all_mean_table_nonsig(z,r) = mean(all_non_sig_struct(z).data,'omitnan');
            try all_pval(z,r) = ranksum(seg_EVE_roi(r).data,all_non_sig_struct(z).data);
            catch error
                all_pval(z,r) = NaN;
            end
            num_non_selected(z,r) = length(ns_seg_EVE_roi(r).data);
            sig_slopes{z,r} = seg_EVE_roi(r).data;
            nonsig_slopes{z,r} = all_non_sig_struct(z).data;
        end
        all_pval_gm(z) = ranksum(all_vals_sig(z).seg,all_non_sig_struct(z).data);
        all_pval_wm(z) = ranksum(all_vals_sig(z).wm,all_non_sig_struct(z).data);
            
        % scatter plot figure for comparing nonselected features to
        % selected features in each anatomical region
        figure(z);clf
        set(gcf,'Color','white','Position',[723   523   577   320])
        hold on
        qqq= 0;
        a1 = scatter(ones(1,length(all_non_sig_struct(z).data)),all_non_sig_struct(z).data,10,'MarkerEdgeColor',color1,'MarkerFaceColor',color1,'jitter','on')
        plot([0.75 1.25],[nanmedian(all_non_sig_struct(z).data) nanmedian(all_non_sig_struct(z).data)],'k-','LineWidth',1)
        for r = 1:8 %regions
            if r<4
                which_color = gm_color_1;
            elseif r<9
                which_color = wm_color_2;
            else
                which_color = color1;
            end
            try [p99(z).data(r) p99(z).h(r) p99(z).stats(r)] = ranksum(seg_EVE_roi(r).data,all_non_sig_struct(z).data);
                    num_selected(z,r) = length(seg_EVE_roi(r).data);
                    
                    qqq = r+2;
                    scatter((r+2)*ones(1,length(seg_EVE_roi(r).data)),seg_EVE_roi(r).data,10,'MarkerEdgeColor',which_color,'MarkerFaceColor',which_color,'jitter','on')
                    plot([r+1.75 r+2.25],[nanmedian(seg_EVE_roi(r).data) nanmedian(seg_EVE_roi(r).data)],'k-','LineWidth',1)

            catch ME
                p99(z).data(r) = 1;
            end
        end
        for r = 1:11
            try max_regions(r,1) = nanmax(seg_EVE_roi(r).data);
            catch ME
                max_regions(r,1) = 0;
            end
        end
        y_val = nanmax([max_regions;all_non_sig_struct(z).data(:)]);
        num_sig = sum(p99(z).data<(0.05/64));
        which_sig = find(p99(z).data<(0.05/64));
        for j = 1:num_sig
            stop_sig = which_sig(j)+2;
            plot([1,stop_sig], [y_val*1.1.^j,y_val*1.1.^j], '-k', 'LineWidth',.5)
            if p99(z).data(stop_sig-2)<(0.001./64)
            plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), '+k','MarkerSize',4)
            elseif p99(z).data(stop_sig-2)<(0.01./64)
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), 'dk','MarkerSize',4)
            elseif p99(z).data(stop_sig-2)<(0.05./64)
                plot(mean([1,stop_sig]), y_val*(1.1.^j+0.05), '*k','MarkerSize',4)
            end
        end
        hold off
        xticks([1,3:qqq])
        xticklabels(['Non-selected',plot_x_label])
        title(sprintf('Metric = %s Frequency = %s lobar ROI',qexp_vs_pow{metric},freq_bands{freq}))
    end
    
    if strcmp(qexp_vs_pow{metric},'qexp')
        sig_ch_all_Q = sig_ch_all;
        sig_ch_names_all_Q = sig_ch_names_all;
    else
        sig_ch_all_P = sig_ch_all;
        sig_ch_names_all_P = sig_ch_names_all;
    end
end

%% 5. anatomical region heat maps (MH)
% create heat map summary of section 4
close all
addAsterisk = 1;

load color_bar.mat
fband_labs = {'{\theta}/{\alpha}','\beta','L_{\gamma}','H_{\gamma}'};
plot_x_label = {'Frontal','Tmp-parietal','Paralimbic',...
                'Thal.-cortical','Frontal Assoc.','Tmp.Assoc.',...
                'Paralimbic','Commissural'};
color_gm = '0.02,0.62,0.584';
color_wm = '0.439,0.024,0.631';
region_labs = {strcat('\color[rgb]{',color_gm,'}Frontal'),...
    strcat('\color[rgb]{',color_gm,'}Tmp-parietal'),...
    strcat('\color[rgb]{',color_gm,'}Paralimbic'),...
    strcat('\color[rgb]{',color_wm,'}Thal-cortical'),...
    strcat('\color[rgb]{',color_wm,'}Frontal Assoc.'),...
    strcat('\color[rgb]{',color_wm,'}Tmp Assoc.'),...
    strcat('\color[rgb]{',color_wm,'}Paralimbic'),...
    strcat('\color[rgb]{',color_wm,'}Commissural')};
sigStars = ["$\ast$", "$\ast\ast$", "$\ast\ast\ast$"];

sig_all_median_table = NaN*zeros(size(all_mean_table)); % 4 fbands * 2 metrics x 8 regions
for z = 1:size(sig_all_median_table,1) % 4 fbands * 2 metrics 
    % which node-features have distributions significantly different from non-significant features
    num_sig(z) = sum(p99(z).data<(0.05/64));
    which_sig = find(p99(z).data<(0.05/64));
    sig_all_median_table(z,which_sig) = all_median_table(z,which_sig);
end

if save_heatmap_data
    %save data for figures w diff thresholds
    t = num2str(percThresh);
    save([savedir savefname],"sig_all_median_table","all_median_table","all_mean_table","all_median_table_nonsig","all_mean_table_nonsig","p99")
end

figure(3); clf;
set(gcf,'Color','white','Units','inches','Position',[5,5,15.5,3.3]);
figpos = get(gcf, 'Position');
set(gcf,'PaperOrientation','landscape','PaperUnits', 'inches','PaperSize', flip(figpos(3:4)));
sgtitle(['Threshold = ' num2str(percThresh)])
a1 = subplot(1,2,1);
h1 = imagesc(sig_all_median_table(1:4,:));
set(h1, 'AlphaData',~isnan(sig_all_median_table(1:4,:)))
a1.XAxis.FontSize = 7;
set(a1,'FontSize',9)
title('PLV average slopes')
yticks([1:length(fband_labs)])
xticks([1:length(region_labs)])
yticklabels(fband_labs)
xticklabels(region_labs)
xtickangle(-30)
colorbar;
colormap(flip(color_bar));
clim([-200,200])

a2 = subplot(1,2,2);
h2 = imagesc(sig_all_median_table(5:8,:));
set(h2, 'AlphaData',~isnan(sig_all_median_table(5:8,:)))
a2.XAxis.FontSize = 7;
set(a2,'FontSize',9)
title('Power average slopes')
yticks([1:length(fband_labs)])
xticks([1:length(region_labs)])
yticklabels(fband_labs)
xticklabels(region_labs)
colorbar;
colormap(flip(color_bar));
xtickangle(-30)
clim([-200,200])

if addAsterisk
    % sig stars
    for z = 1:8 
        which_sig = find(p99(z).data<(0.05/64));
        if z <=4
            f = z;
            subplotNum = 1;
        else
            f = z-4;
            subplotNum = 2;
        end
    
        for j = 1:num_sig(z)
            disp(['pval: ' num2str(p99(z).data(which_sig(j)))])
            figure(3);
            subplot(1,2,subplotNum)
            if p99(z).data(which_sig(j))<(0.001./64)
                text(which_sig(j),f, sigStars(3), 'Interpreter', 'latex','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            elseif p99(z).data(which_sig(j))<(0.01./64)
                text(which_sig(j),f, sigStars(2), 'Interpreter', 'latex','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            elseif p99(z).data(which_sig(j))<(0.05./64)
                text(which_sig(j),f, sigStars(1), 'Interpreter', 'latex','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');
            end
        end
    end
end

%% 6. Intrinsic communicability of anatomical regions and subregions
% plot distribution of z-scored features values in each region across all trials

% run this section for each composite region of interest
regIdx = 1; % 1 = FrGM 2 = TPGM 3 = PLGM 4 = TCWM 5 = FAWM 6 = TAWM 7 = PLWM 8 = ComWM

%load patient loc data & fix issue with subj 18
load patient_loc_120623.mat
goodchs = patient_loc(1).session(18).type~=0;
patient_loc(1).session(18).names = patient_loc(1).session(18).names(goodchs);
patient_loc(1).session(18).gm_wm_rois = patient_loc(1).session(18).gm_wm_rois(goodchs);
patient_loc(1).session(18).gm_wm_labels = patient_loc(1).session(18).gm_wm_labels(goodchs);
patient_loc(1).session(18).seg_ID = patient_loc(1).session(18).seg_ID(goodchs);

aQv = avgQ_allFeats;
fIDq = fIDq_allFeats;
gchQ = gchQ_allFeats;
ssIDq = ssIDq_allFeats;

code_regions_qexp = zeros(size(aQv));
%generate code_regions variable
for i = 1:length(avgQ_allFeats)
    iSubj = ssIDq(i);
    chIdx = find(strcmp(patient_loc(1).session(iSubj).names, gchQ(i)));
    if ~isempty(chIdx)
        code_regions_qexp(i) = patient_loc(1).session(iSubj).gm_wm_rois(chIdx); 
    end
end

region_labels = {'FrGM','TPGM','PLGM','TCWM','FAWM','TAWM','PLWM','ComWM'};

% plot subregions in different colors - get subregion labels
% gm & wm ROIs
EVE_roi_WM_csv = EVE_roi_WM_csv{:,2};
EVE_roi_WM_csv{end+1} = 'n/a';
seg_ID_list = [frontal_cortex,temporoparietal_cortex,paralimbic_GM,'n/a'];

%subregion codes
FrGM_sub = [1:7]';
TPGM_sub = [8:15]';
PLGM_sub = [16:22]';
TCWM_sub = [36,96; 7,69; 6,68; 38,98; 37,97; 35,95];
FAWM_sub = [46,106; 5,67; 21,82; 4,66; 24,85; 3,65; 43,103; 42,102];
TAWM_sub = [8,70; 12,74; 19,80; 44,104; 20,81; 18,79; 1,64; 45,105];
PLWM_sub = [39,99; 40,100; 41,101; 11,73; 47,107];
ComWM_sub = [52,112; 51,111; 53,113];
subregion_codes = {FrGM_sub, TPGM_sub, PLGM_sub, TCWM_sub, FAWM_sub, TAWM_sub,PLWM_sub,ComWM_sub};

% label
subregionIdx = zeros(size(gchQ));
gchQ_gmwmlabel = cell(size(gchQ));
for i = 1:length(gchQ)
    iSubj = ssIDq(i);
    iReg = code_regions_qexp(i);
    chIdx = find(strcmp(patient_loc(1).session(iSubj).names, gchQ(i)));
    if ~isempty(chIdx)
    if iReg<=3 %gm region
        sub_code = patient_loc(1).session(iSubj).seg_ID(chIdx);
        if sub_code>0
            subregionIdx(i) = find(subregion_codes{iReg}==sub_code);
        end
    else %wm region
        sub_code = patient_loc(1).session(iSubj).eve(chIdx);
        if sub_code>0
            [subi,subj] = find(subregion_codes{iReg}==sub_code);
            if isempty(subi)
                x=1;
            else
                subregionIdx(i) = subi;
            end
        end
    end
    end
end

% set up colors for plot
color_assignment = [1,0,0; .3,.73,.05; .05,.4,.73; .9,.53,.1; 1,0,1; .27,.8,.85; .8,0,.5; 1,1,0]; %white background
xtick_vals = [1,2.2,3.4,4.6,5.8,7,8.2,9.4,10.6,11.8];
subregion_labels = {'All Subregions','1','2','3','4','5','6','7','8'};
color_vect = zeros(length(subregionIdx),3);
for i = 1:length(subregionIdx)
    if subregionIdx(i)>0
        color_vect(i,:)=color_assignment(subregionIdx(i),:);
    end
end

%plot intrinsic qexp for subregions
plot_feats = aQv(code_regions_qexp == regIdx);
plot_subregionIdx = subregionIdx(code_regions_qexp == regIdx);
plot_color_vect = color_vect(code_regions_qexp == regIdx,:);
plot_feats(plot_subregionIdx==0) = [];
plot_color_vect(plot_subregionIdx==0,:) = [];
plot_subregionIdx(plot_subregionIdx==0) = [];

figure(3)
hold on
set(gcf,"Color","white")
numSubs = height(subregion_codes{regIdx});
for i = 1:numSubs+1 %subregion
    figure(3)
    if i==1
        y = plot_feats;
        x = ones(length(y),1);
        swarmchart(x*xtick_vals(i),y,x*20,plot_color_vect,'o',"MarkerFaceColor","flat")
        box off
    else
        y = plot_feats(plot_subregionIdx==i-1);
        x = ones(length(y),1);
        swarmchart(x*xtick_vals(i),y,x*20,color_assignment(i-1,:),'o',"MarkerFaceColor","flat")
        box off
    end
    ylim([-1.5,2.5])
    xlim([0,xtick_vals(i)+1])
    xticks(xtick_vals)
    xticklabels(subregion_labels)
    if blackBackground
        yline(0,'w','LineStyle','--',"LineWidth",1.5)
    else
        yline(0,'k','LineStyle','--',"LineWidth",1.5)
    end
    if ~isempty(y)
        p(1,i) = signrank(y,0);
        if p(1,i) < 0.001
            scatter(xtick_vals(i),2.25,'*','k')
        elseif p(1,i) < 0.01
            scatter(xtick_vals(i),2.25,'+','k')
        elseif p(1,i) < 0.05
            scatter(xtick_vals(i),2.25,'o','k')
        end
    end
end
