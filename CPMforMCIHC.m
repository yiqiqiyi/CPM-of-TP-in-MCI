% This code is modified based on the code of 2015 Xilin Shen and Emily Finn.

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,
% Papademetris X & Constable RT. (2015). Functional connectome
% fingerprinting: Identifying individuals using patterns of brain
% connectivity. Nature Neuroscience 18, 1664-1671.

% This code provides a framework for implementing functional
% connectivity-based behavioral prediction in a leave-one-subject-out
% cross-validation scheme, as described in Finn, Shen et al 2015 (see above
% for full reference). The first input ('all_mats') is a pre-calculated
% MxMxN matrix containing all individual-subject connectivity matrices,
% where M = number of nodes in the chosen brain atlas and N = number of
% subjects. Each element (i,j,k) in these matrices represents the
% correlation between the BOLD timecourses of nodes i and j in subject k
% during a single fMRI session. The second input ('all_behav') is the
% Nx1 vector of scores for the behavior of interest for all subjects.

% As in the reference paper, the predictive power of the model is assessed
% via correlation between predicted and observed scores across all
% subjects. Note that this assumes normal or near-normal distributions for
% both vectors, and does not assess absolute accuracy of predictions (only
% relative accuracy within the sample). It is recommended to explore
% additional/alternative metrics for assessing predictive power, such as
% prediction error sum of squares or prediction r^2.

clc;clear;
%exclude data based on headmotion
x1=importdata('HeadMotion1.tsv');
HeadMotion1=load('HeadMotion1.mat') %load headmation
x2=importdata('HeadMotion2.tsv');
HeadMotion2=load('HeadMotion2.mat') %load headmation
x.data=[x1.data;x2.data];
x.textdata=[x1.textdata(:,1);x2.textdata(2:end,1)];
HeadMotion=[HeadMotion1.HeadMotion;HeadMotion2.HeadMotion];
headnum=x.textdata(2:end);  %get the rank of headmotion
clear headidx
%rename the individual
for hn=1:size(headnum,1)
    hmi=headnum{hn};
    if length(hmi)==3
    headidx(hn,1)=str2num(['10',hmi(2:end)]);
    elseif length(hmi)==4
    headidx(hn,1)=str2num(['1',hmi(2:end)]);
    else
    headidx(hn,1)=str2num(hmi(2:end));
    end
end

meanJK=find(HeadMotion(:,20)>0.25);
maxheadmotion=3;
maxhemo=[];
maxhemo=find(HeadMotion(:,1)>maxheadmotion|HeadMotion(:,2)>maxheadmotion|HeadMotion(:,3)>maxheadmotion|...
    HeadMotion(:,4)>maxheadmotion|HeadMotion(:,5)>maxheadmotion|HeadMotion(:,6)>maxheadmotion);
subj_excl_HM=union(meanJK,maxhemo);
subj_excl=headidx(subj_excl_HM);

load CPM_gloCWF2_1.mat；
load CPM_gloCWF2_2.mat；
rest_1_mats=cat(3,rest_1mats_1,rest_1mats_2)；

trust_tet=readtable('trust6_CPM.xlsx','ReadVariableNames',1)
PMAT_CR=table2array(trust_tet(:,3));

%% ------------ INPUTS -------------------
%Subject str2num
clear sub_name
for s=1:size(trust_tet,1)
    sub_num=num2str(table2array(trust_tet(s,1)));
    if sub_num(1)=='1'
    if sub_num(2)=='0'
        sub_name(s)=str2num(['10',sub_num(3:4)]);
    else
        sub_name(s)=str2num(['1',sub_num(2:4)]);
    end
    else
        sub_name(s)=str2num([sub_num(1:4)]);
    end
end
%select group
%MCI
mci_grp=find(table2array(trust_tet(:,11))==2);
%HC
hc_grp= find(table2array(trust_tet(:,11))==1);
%excluded the subject with excessive headmotion or sleep in scanning or not
%not matching the rule of one-shot TG
exclu_idx=intersect(sub_name,subj_excl);
exclu_idxnum=find(ismember(sub_name,exclu_idx)==1);
remain_subj=sub_name;
remain_subj(exclu_idxnum)=[];
manual_del=[8 106 114 116 147 169 175];
mci_allgrp1=setdiff(mci_grp,[exclu_idxnum,manual_del]);
hc_allgrp1=setdiff(hc_grp,[exclu_idxnum,manual_del]);
%select group
allgrp=[hc_allgrp1]
%select MCI data
mci_rest_1_mats=rest_1_mats(:,:,allgrp);
mci_behavior=PMAT_CR(allgrp);
mci_trust_tet=trust_tet(allgrp,:);
%%
%remove cerreblum
load /usrdir1/user02/CBM/basedon_160altas/net2name.mat
lobenum=net2name;
    all_mats  = mci_rest_1_mats;
    relobe=[];
     relobe=find(lobenum(:,1)==6);

    renode=relobe;
    all_mats(renode,:,:)=0;
    all_mats(:,renode,:)=0;

%input
all_behav = mci_behavior;
allage=[];
for a=1:length(all_behav)
allage(a,1)=str2num(cell2mat(table2array(trust_tet(allgrp(a),6))));
end
alledu=table2array(trust_tet(allgrp,8));
allsex=table2array(trust_tet(allgrp,7));
clear allhead
remain_subj_head=allgrp;
% remain_subj_head(exclu_idxnum1)=[];
for p=1:length(remain_subj_head)
headidx_trust=find(headidx==table2array(trust_tet(remain_subj_head(p),1)));
allhead(p,1)=HeadMotion(headidx_trust,20);
end

% threshold for feature selection
thresh = 0.01;

% ---------------------------------------

no_sub = size(all_mats,3);
no_node = size(all_mats,1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);
behav_pred = zeros(no_sub,1);

for leftout = 1:no_sub;
    fprintf('\n Leaving out subj # %6.3f',leftout);
    
    % leave out subject from matrices and behavior
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    
    train_behav = all_behav;
    train_behav(leftout) = [];
    
    % correlate all edges with behavior using partial correlation
    train_age=allage;
    train_age(leftout) =[];
    train_edu=alledu;
    train_edu(leftout) =[];
    train_sex=allsex;
    train_sex(leftout) =[];
    train_head=allhead;
    train_head(leftout) =[];
    %normalize data
     [ztrain_age,mtrain_age,sitrain_age]=zscore(train_age);
     [ztrain_edu,mtrain_edu,sitrain_edu]=zscore(train_edu);
      [ztrain_sex,mtrain_sex,sitrain_sex]=zscore(train_sex);
     [ztrain_head,mtrain_head,sitrain_head]=zscore(train_head);
      [ztrain_behav,mtrain_behav,sitrain_behav]=zscore(train_behav);
      [ztrain_vcts,mtrain_vcts,sitrain_vcts]=zscore(train_vcts');
%        
%     % correlate all edges with behavior using rank correlation
     [r_mat, p_mat] = partialcorr(ztrain_vcts, ztrain_behav, [ztrain_age,ztrain_edu,ztrain_sex,ztrain_head], 'type', 'Spearman');
    ztrain_mats=reshape(ztrain_vcts,[],no_node,no_node);
        
    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    
    % set threshold and define masks 
    pos_mask = zeros(no_node, no_node);
    neg_mask = zeros(no_node, no_node);
    
    
    pos_edge = find( r_mat >0 & p_mat < thresh);
    neg_edge = find( r_mat <0 & p_mat < thresh);
    
    pos_mask(pos_edge) = 1;
    neg_mask(neg_edge) = 1;
        
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = zeros(no_sub-1,1);
    
    for ss = 1:size(train_sumpos);
        
        train_sumpos(ss) = sum(nansum(squeeze(ztrain_mats(ss,:,:)).*pos_mask))/2;
        train_sumneg(ss) = sum(nansum(squeeze(ztrain_mats(ss,:,:)).*neg_mask))/2;
    end
    
    % build model on TRAIN subs
    % combining both postive and negative features
    
    
    b = regress(ztrain_behav, [train_sumpos, train_sumneg, ones(no_sub-1,1)]);
    fit_pos = polyfit(train_sumpos, ztrain_behav,1);
    fit_neg = polyfit(train_sumneg, ztrain_behav,1);

    % run model on TEST sub
    
    test_mat = all_mats(:,:,leftout);
    test_vcts = reshape(test_mat,[],size(test_mat,3));
    for zn=1:length(test_vcts)
       ztest_vcts(zn,1) =(test_vcts(zn)-mtrain_vcts(zn))./sitrain_vcts(zn);
    end
    ztest_mat=reshape(ztest_vcts,no_node,no_node);
    test_sumpos = sum(nansum(ztest_mat.*pos_mask))/2;
    test_sumneg = sum(nansum(ztest_mat.*neg_mask))/2;
    
    behav_pred(leftout) = b(1)*test_sumpos + b(2)*test_sumneg + b(3);
    behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    leav_pos(:,:,leftout)=pos_mask;
    leav_neg(:,:,leftout)=neg_mask;
    testnormal(leftout,1)=(all_behav(leftout)-mtrain_behav)/sitrain_behav;
    corpre(leftout,1)=corr(train_behav,train_sumpos);
    corpre(leftout,2)=corr(train_behav,train_sumneg);
end

% Estimating R
[R_pos, P_pos] = corr(behav_pred_pos,testnormal, 'type', 'Spearman')
[R_neg, P_neg] = corr(behav_pred_neg,testnormal, 'type', 'Spearman')
[R_com, P_com] = corr(behav_pred,testnormal, 'type', 'Spearman')
% Estimating MSE
MSE_pos=sum((behav_pred_pos-testnormal).^2)/(length(testnormal)-length(fit_pos)-1)
MSE_neg=sum((behav_pred_neg-testnormal).^2)/(length(testnormal)-length(fit_neg)-1)
MSE_sum=sum((behav_pred-testnormal).^2)/(length(testnormal)-length(b)-1)
% figure
figure(1); plot(testnormal,behav_pred_pos,'r.'); lsline;title(['Positive model (r = ',num2str(round(R_pos,2)),')'])
figure(2); plot(testnormal,behav_pred_neg,'b.'); lsline;title(['Negative model (r = ',num2str(round(R_neg,2)),')'])
figure(3); plot(testnormal,behav_pred,'k.'); lsline;title(['Com model (r = ',num2str(round(R_com,2)),')'])
