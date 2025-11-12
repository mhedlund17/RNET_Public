% make time-resolved PLV figures with significance tests
% averages ACR-to-all and ACR-to-FrGM and ACR-to-TPGM pairs over all pairs


%% load data & settings
% Set composite region for each region of interest (integer 1-8)
    % 1 = Frontal GM
    % 2 = Temporoparietal GM
    % 3 = Paralimbic GM
    % 4 = Thalamocortical WM
    % 5 = Frontal Association WM 
    % 6 = Temporal Association WM
    % 7 = Paralimbic WM
    % 8 = Commissural WM
reg1_composite = 4;
reg2_composite = 1;

% Set anatomical or functional subregion of each region (empty or character array)
% leave empty (=[]) to define region as composite region
    % 'ACR' = Anterior Corona Radiata in TCWM
    % 'SCR' = Superior Corona Radiata in TCWM
    % 'PTR' = Posterior Thalamic Radiations in TCWM
    % 'PrWM' = Precentral WM in TCWM
    % 'SFG' = Suprior Frontal Gyrus in FrGM
    % 'MFG' = Middle Frontal Gyrus in FrGM
    % 'IFG' = Inferior Frontal Gyrus in FrGM
    % 'PrG' = Precentral Gyrus in FrGM
    % 'SMA' = Supplementary Motor Area in FrGM (functional)
    % 'PMA' = PreMotor Area in FrGM (functional)
    % 'dmPFC' = Dorsomedial Prefrontal Cortex in FrGM (functional)
    % 'dlPFC' = Dorsolateral Prefrontal Cortex in FrGM (functional)
reg1_subreg = ['ACR'];
reg2_subreg = ['MFG'];

% regPref = 'ACR'; %'ACR' 'SCR' 'Precent' 'PTR'
conditions = [1, 2, 3]; %1 = fast trials, 2 = slow trials, 3 = all trials
Ncond = length(conditions);
% fig_odir = ''; % directory to save figures in
stsubj = 1; %start subject

db = CCDTdatabase;
Nsubj = height(db);
ddir = 'test_output\'; % directory to time-resolved PLV data 
load patient_loc_MH_071025.mat %anatomical info
fbandlabels = {'Theta/Alpha', 'Beta', 'Low Gamma', 'High Gamma'};
plotSingSubj = 0;
plotMultiSubj = 1;

% gm seg ROI
frontal_cortex = {'SFG', 'MFG', 'IFG', 'precentral', 'GRe', 'MOrG', 'POrG'};
temporoparietal_cortex = {'MTG', 'ITG', 'STG', 'FuG', 'AnG', 'SMG', 'SPL','postcentral'};
paralimbic_GM = {'Hippocampus','Amygdala','PHG','cingulate','LiG','Ent','insula'};
seg_ID_list = [frontal_cortex,temporoparietal_cortex,paralimbic_GM,'n/a'];

%% get pairs between specified regions
%initiate multisubj channel pair structure
AllSubjData.pairs = cell(Nsubj,1);

for isubj = stsubj:Nsubj
    disp(num2str(isubj))

    % get pairs
    Nch = height(patient_loc(1).session(isubj).names);
    cgch = 1:Nch; % number of good channels (from loadCCDTdata)
    chprs = nchoosek(cgch,2); % all pairs
    Npr = size(chprs,1);
    inBrain = patient_loc(1).session(isubj).type~=0; %channels labeled inside brain - fixes issue with subj 18 patient_loc file.

    for iReg = 1:2
        subregChs = [];
        if iReg ==1
            compositeLabel = reg1_composite;
            subregLabel = reg1_subreg;
        else
            compositeLabel = reg2_composite;
            subregLabel = reg2_subreg;
        end
        
        % find composite region channels
        compositeChs = patient_loc(1).session(isubj).gm_wm_rois(inBrain)==compositeLabel;
        
        %find subregion channels
        if isempty(subregLabel) % region is composite region
            subregChs = compositeChs;
        else
            % check that composite channel matches subregion
            if (strcmp(subregLabel,'ACR')||strcmp(subregLabel,'SCR')||strcmp(subregLabel,'PTR')||...
                strcmp(subregLabel,'PrWM')) && compositeLabel ~= 4
                disp(['Composite region for subregion ' num2str(iReg) ' does not match. ' ...
                    subregLabel ' is in TCWM. Change the composite region to 4.'])
                break
            end
            if (strcmp(subregLabel,'SFG')||strcmp(subregLabel,'MFG')||strcmp(subregLabel,'IFG')||...
                strcmp(subregLabel,'PrG')||strcmp(subregLabel,'SMA')||strcmp(subregLabel,'PMA')||...
                strcmp(subregLabel,'dmPFC')||strcmp(subregLabel,'dlPFC')) && (compositeLabel ~= 1)
                disp(['Composite region for subregion ' num2str(iReg) ' does not match. \n' ...
                    subregLabel ' is in FrGM. Change the composite region to 1.'])
                break
            end

            % find subregion channels
            if strcmp(subregLabel,'SFG')
                subregChs = compositeChs & patient_loc(1).session(isubj).seg_ID(inBrain)==1; %SFG seg_ID_idx(1)
            elseif strcmp(subregLabel,'MFG')
                subregChs = compositeChs & patient_loc(1).session(isubj).seg_ID(inBrain)==2; %MFG seg_ID_idx(2)
            elseif strcmp(subregLabel,'IFG')
                subregChs = compositeChs & patient_loc(1).session(isubj).seg_ID(inBrain)==3; %IFG seg_ID_idx(3)
            elseif strcmp(subregLabel,'PrG')
                subregChs = compositeChs & patient_loc(1).session(isubj).seg_ID(inBrain)==4; %PrG seg_ID_idx(4)
            elseif strcmp(subregLabel,'SMA')
                subregChs = patient_loc(1).session(isubj).inSMA(inBrain);
            elseif strcmp(subregLabel,'PMA')
                subregChs = patient_loc(1).session(isubj).inPMA(inBrain);
            elseif strcmp(subregLabel,'dmPFC')
                subregChs = patient_loc(1).session(isubj).indmPFC(inBrain);
            elseif strcmp(subregLabel,'dlPFC')
                subregChs = patient_loc(1).session(isubj).indlPFC(inBrain);
            elseif strcmp(subregLabel,'ACR')
                subregChs = compositeChs & startsWith(patient_loc(1).session(isubj).wm_label,'Anterior_corona_radiata');
            elseif strcmp(subregLabel,'SCR')
                subregChs = compositeChs & startsWith(patient_loc(1).session(isubj).wm_label,'Superior_corona_radiata');
            elseif strcmp(subregLabel,'PrWM')
                subregChs = compositeChs & startsWith(patient_loc(1).session(isubj).wm_label,'PRECENTRAL__WM');
            elseif strcmp(subregLabel,'PTR')
                subregChs = compositeChs & startsWith(patient_loc(1).session(isubj).wm_label,'Posterior_thalamic_radiation'); 
            end
        end
        % get channel index
        chsIdx = find(subregChs);
        if iReg==1
            chsIdx_1 = chsIdx;
        else
            chsIdx_2 = chsIdx;
        end
    end

    % get all pairs between both regions
    allPairsIdx = [];
    for i = 1:length(chsIdx_1)
        for j = 1:length(chsIdx_2)
            prIdx = find((chprs(:,1)==chsIdx_1(i) & chprs(:,2)==chsIdx_2(j)) | (chprs(:,2)==chsIdx_1(i) & chprs(:,1)==chsIdx_2(j)));
            allPairsIdx = [allPairsIdx;prIdx]; %concatenate all pairs between regions
        end
    end
    AllSubjData.pairs{isubj} = chprs(unique(allPairsIdx),:);
    AllSubjData.pairIdx{isubj} = allPairsIdx;   
end
%% Plot time-resolved PLV differences between fast and slow conditions
%initiate multisubj structure for pairs structure
for f = 1:4
    AllSubjData.fbands(f).pairData = cell(Ncond,Nsubj); %average across all chs
end
AllSubjData.t = cell(Nsubj,1); %time variable

disp('Getting Single Subject Data')
flist = dir(ddir);
for isubj = stsubj:Nsubj %get pairwise data for all subjects
    disp(isubj)
    if ~isempty(AllSubjData.pairs{isubj}) %skip subjs with no pairs
    %load trPLV data
    %edit this according to the file naming system. Below code assumes
    %file begins with subject and session name.
    fname_prefix = [db{isubj,1} '_' db{isubj,2}];
    for i = 1:height(flist)
        if startsWith(flist(i).name,fname_prefix)
            fname = flist(i).name;
        end
    end
    load([ddir fname]);
    
    Nch = TRstruct.Nch; Nsamp = TRstruct.Nsamp; 
    Nfbands = length(TRstruct.fbands); win = TRstruct.win;
    cgch = 1:Nch; % number of good channels (from loadCCDTdata)
    chprs = nchoosek(cgch,2); % all pairs
    AllSubjData.t{isubj} = TRstruct.t; 

    %get idx of each pair --> idx are same as AllSubjData.pairIdx, but in ascending order
    prIdx = [];
    for i = 1:height(AllSubjData.pairs{isubj}) %get acr-to-frgm pair idx
        prIdx(i) = find(chprs(:,1)==AllSubjData.pairs{isubj}(i,1) & chprs(:,2)==AllSubjData.pairs{isubj}(i,2));
    end
    
    %get trplv data for pairs
    baseWin = [-750 -500]; %prep/antic period window is [-500,0]ms, so baseline is 250ms window before that window
    baseWinIdx = find(TRstruct.t>baseWin(1) & TRstruct.t<baseWin(2));
    for f = 1:Nfbands
        for ll = 1:Ncond
            c = conditions(ll);
            AllSubjData.fbands(f).pairData{c,isubj} = TRstruct.fbands(f).atplv(prIdx,:,c);
        end
        % normalize to common baseline - avg of fast and slow conditions. condition 1 is fast, condition 2 is slow
        base_mu_fast = mean(AllSubjData.fbands(f).pairData{1,isubj}(:,baseWinIdx),2);
        base_mu_slow = mean(AllSubjData.fbands(f).pairData{2,isubj}(:,baseWinIdx),2);
        base_mu = (base_mu_slow+base_mu_fast)./2;
        base_sig_fast = std(AllSubjData.fbands(f).pairData{1,isubj}(:,baseWinIdx),0,2);
        base_sig_slow = std(AllSubjData.fbands(f).pairData{2,isubj}(:,baseWinIdx),0,2);
        base_sig = (base_sig_slow+base_sig_fast)./2;            
        c=1; AllSubjData.fbands(f).pairData{c,isubj} = (AllSubjData.fbands(f).pairData{c,isubj}-base_mu)./base_sig;  
        c=2; AllSubjData.fbands(f).pairData{c,isubj} = (AllSubjData.fbands(f).pairData{c,isubj}-base_mu)./base_sig;  
    end
    clear TRstruct %remove from memory

    % plot individual subject, fast and slow conditions
    if plotSingSubj
        figure;
        set(gcf, 'Position', get(0, 'Screensize'),'Color','white');
        sgtitle(fname_prefix,'Interpreter','none')
        for f = 1:Nfbands
            subplot(2,2,f)
            hold on
            for ll = 1:2 % fast and slow conditions
                c = conditions(ll);
                plot(AllSubjData.t{isubj},mean(AllSubjData.fbands(f).pairData{c,isubj},1)) %average trace
            end
            ylim([-3,3])
            xlim([win(1)+100,win(2)-100]) %100ms buffer
            xlabel('Time (ms)')
            ylabel('Average Normalized PLV')
            title([fbandlabels{f}])
            xline(0,'k',LineWidth=2,Alpha=.5)
            xline(-500,'k--',LineWidth=2,Alpha=.5)
            legend('fast','slow')
        end
    end
    
    if length(AllSubjData.t{isubj})<2000
        %some subjs have different sampling rate. 
        % set t to shorter sampling rate - will downsample other subjects in multisubj step
        t = AllSubjData.t{isubj}; 
    end
    close all
    end
end %subj loop

% plot pairs across all subjects
if plotMultiSubj
    %combine pairs from all subjects & downsample subjs w fs=1024
    cumulative_pairData = cell(Nfbands,Ncond);
    for isubj = 1:Nsubj
    if ~isempty(AllSubjData.pairs{isubj}) %skip subjs with no pairs
        for f = 1:Nfbands
            for ll = 1:Ncond
                c = conditions(ll);
                if width(AllSubjData.fbands(f).pairData{c,isubj}) > length(t)
                    %downsample 1024 sampling rate to 512
                    AllSubjData.fbands(f).pairData{c,isubj} = downsample(AllSubjData.fbands(f).pairData{c,isubj}',2)';
                end
                %concatenate data from all subjs
                cumulative_pairData{f,c} = vertcat(cumulative_pairData{f,c},AllSubjData.fbands(f).pairData{c,isubj});
            end
        end
    end
    end %subj

% fast vs slow statistical test
yValue= 0.7; yValue2=0.65; sample_min = 20; 
figvar = figure();
% figpos = [10.5895    6.1158    9.4947    5.1368];
% set(figvar,'Units','inches'); set(figvar,'Position',figpos);
set(figvar, 'Position', get(0, 'Screensize'),'Color','white');
set(figvar,'Color','white','Units','inches','PaperUnits','inches');
% set(figvar,'PaperSize',figpos(3:4),'Renderer','painters');
sgtitle('PLV Across Subjects')
for f = 1:Nfbands
    %find significant difference between fast and slow conditions at each sample in trPLV analysis window
    num_comp = length(find(t<500 & t>-750))*4; %correct for number of samples, x 4 fbands
    corr_alpha = .05/num_comp; %bonferoni multiple comparison correction 
    p_signrank = zeros(size(t)); h_signrank = zeros(size(t));
    for isamp = 1:length(t)
        [p_signrank(isamp), h_signrank(isamp)] = signrank(cumulative_pairData{f,1}(:,isamp),cumulative_pairData{f,2}(:,isamp),'alpha',corr_alpha);
    end

    subplot(2,2,f)
    hold on
    % Plot significant p-values as a horizontal line:
    % plot all significant samples
    prevIndex = 0; % To handle the start of a new segment
    for i = 1:length(h_signrank)
        if h_signrank(i)
            % If this is part of a significant segment, plot a horizontal line
            if prevIndex == 0 % start of a new significant segment
                prevIndex = i; % Update the starting index of the segment
            end
        else
            % If we reach a non-significant value and we have a significant segment
            if prevIndex > 0 
                % Plot horizontal line for the segment
                plot(t(prevIndex:i-1), repmat(yValue2, 1, i-prevIndex), 'k','LineWidth',0.75);
                prevIndex = 0; % Reset for the next segment
            end
        end
    end
    % Handle case if the last values were significant
    if prevIndex > 0
        plot(t(prevIndex:end), repmat(yValue2, 1, length(t)-prevIndex+1), 'k','LineWidth',0.75);
    end

    % plot sections of samples with at least sample_min significant samples in a ros
    prevIndex = 0; % To handle the start of a new segment
    for i = 1:length(h_signrank)
        if h_signrank(i)
            % If this is part of a significant segment, plot a horizontal line
            if prevIndex == 0 % start of a new significant segment
                prevIndex = i; % Update the starting index of the segment
            end
        else
            % If we reach a non-significant value and we have a significant segment
            if prevIndex > 0 
                % Plot horizontal line for the segment
                if i - prevIndex >= sample_min
                plot(t(prevIndex:i-1), repmat(yValue, 1, i-prevIndex), 'r','LineWidth',0.75);
                end
                prevIndex = 0; % Reset for the next segment
            end
        end
    end
    % Handle case if the last values were significant
    if prevIndex > 0
        plot(t(prevIndex:end), repmat(yValue, 1, length(t)-prevIndex+1), 'r','LineWidth',0.75);
    end
    
    % plot trPLV
    plot(t,mean(cumulative_pairData{f,1},1),'b','LineWidth',.5) %average trace
    plot(t,mean(cumulative_pairData{f,2},1),'g','LineWidth',.5) %average trace
    ylim([-.5,.75])
    xlim([baseWin(1),500]) %100ms buffer
    xlabel('Time (ms)')
    ylabel('Normalized PLV')
    title([fbandlabels{f} ' Npr = ' num2str(height(cumulative_pairData{f,1}))])
%     title(fbandlabels{f})
    xline(0,'k',LineWidth=0.75,Alpha=.5)
    xline(-500,'k--',LineWidth=0.75,Alpha=.5)
    ax = gca;  ax.XColor = 'k';  ax.YColor = 'k';
    xticks([-500,0]); yticks([-.4,0,.4])
%     fontsize(gca,5,'points')
end
% saveas(gcf,[fig_odir regPref 'plv_allSubjs_meanPrs_signrank_iev' num2str(iev) '.png'])

end %do plot across subjects


%% Connectivity tilt analysis
%store difference in cumulative fast/slow signals. row 1 is theta/alpha band. row 2 is high gamma band.
cumul_fs_diff = cell(2,1); meanDim = 1; % dim = 1xNsamp
cumul_fs_diff{1,1} = mean(cumulative_pairData{1,1}-cumulative_pairData{1,2},meanDim); 
cumul_fs_diff{2,1} = mean(cumulative_pairData{4,1}-cumulative_pairData{4,2},meanDim);

% connectivity tilt plots
plotpoints = 0; %plot each point in boxchart, and draw line connecting between TA and HG boxes
xoffset = 0.15;
sigstary1= 1.65;
sigstary2= -1.5;
num_comp = 1; %number of comparisons. change if analyzing multiple regions - in paper this was 9
alpha = 0.05/num_comp; %multiple comparison correction

figure; 
set(gcf, 'Color','white');
sgtitle('Delta PLV','Interpreter','none')
hold on
yTA = cumul_fs_diff{1,1}; xTA = ones(length(yTA),1)-xoffset;
yHG = cumul_fs_diff{2,1}; xHG = ones(length(yHG),1)+xoffset;
boxchart(xTA,yTA,'BoxWidth',.25,'BoxFaceAlpha',0,'BoxEdgeColor','b','MarkerStyle','o','MarkerSize',3,'MarkerColor','k','LineWidth',0.5)
boxchart(xHG,yHG,'BoxWidth',.25,'BoxFaceAlpha',0,'BoxEdgeColor','r','MarkerStyle','o','MarkerSize',3,'MarkerColor','k','LineWidth',0.5)
yline(0,'--','LineWidth',2,'Alpha',.5)
if plotpoints
    plot([1-xoffset;1+xoffset],[yTA;yHG],'ko','LineStyle','-','Color',[0,0,0,.3],'MarkerSize',3,'MarkerFaceColor','k')  
end
%theta/alpha
[p,h] = signrank(yTA,alpha);
if h
    text(1-xoffset,sigstary1,'*','FontSize',16,'Color','k','HorizontalAlignment','center')
    e = meanEffectSize(yTA,'Effect','cohen');
    if abs(e.Effect) >= 1
    text(1-xoffset,sigstary1-.05,num2str(round(e.Effect,2)),'FontSize',9,'Color','r','HorizontalAlignment','center')
    else 
    text(1-xoffset,sigstary1-.05,num2str(round(e.Effect,2)),'FontSize',9,'Color','k','HorizontalAlignment','center')
    end
end
%high gamma
[p,h] = signrank(yHG,alpha);
if h
    text(1+xoffset,sigstary1,'*','FontSize',16,'Color','k','HorizontalAlignment','center')
    e = meanEffectSize(yHG,'Effect','cohen');
    if abs(e.Effect) >= 1
    text(1+xoffset,sigstary1-.05,num2str(round(e.Effect,2)),'FontSize',9,'Color','r','HorizontalAlignment','center')
    else 
    text(1+xoffset,sigstary1-.05,num2str(round(e.Effect,2)),'FontSize',9,'Color','k','HorizontalAlignment','center')
    end
end
%theta/alpha vs high gamma
[p,h] = signrank(yTA,yHG,alpha);
if h
    text(1,sigstary2,'*','FontSize',16,'Color','k','HorizontalAlignment','center')
    e = meanEffectSize(yTA,yHG,'Effect','cohen');
    if abs(e.Effect) >= 1.8
    text(1,sigstary2-.05,num2str(round(e.Effect,2)),'FontSize',9,'Color','r','HorizontalAlignment','center')
    else 
    text(1,sigstary2-.05,num2str(round(e.Effect,2)),'FontSize',9,'Color','k','HorizontalAlignment','center')
    end
end

legend('Theta/Alpha','High Gamma')
ylabel('Delta PLV')
ylim([-1.75,1.75])
