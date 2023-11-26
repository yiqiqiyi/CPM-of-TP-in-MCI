clc;clear;
% trust_txt=xlsread('//mnt/user02_3/trust5.xlsx');
trust_txt=readtable('//mnt/user02_3/trust5.xlsx','ReadVariableNames',1)
trust_txt([134 138 162],:)=[];
PMAT_CR=trust_txt(:,4);
trust_subj=table2array(trust_txt(:,1));
cd('//mnt/user02_3/resting_analysis/Results/ROISignals_FunImgARglobalCWF//');
for i=1:length(trust_subj)
    sub_num=num2str(trust_subj(i));
    if sub_num(1)=='1'
    if sub_num(2)=='0'
        sub_name=['S',sub_num(3:4)];
    else
        sub_name=['S',sub_num(2:4)];
    end
    else
        sub_name=['S',sub_num(1:4)];
    end
    name{i,1}=sub_name;
    cormatrix=['ROICorrelation_FisherZ_',sub_name,'.mat'];
    load(cormatrix);
    rest_1_mats(:,:,i)=ROICorrelation_FisherZ;
end
 
%  mkdir('/mnt/user02_3/trust_MCI_264')
 cd('//mnt/user02_3/trust_MCI_268');
 save('CPM_gloCWF.mat','PMAT_CR','rest_1_mats');

 
 
 clc;clear;
% trust_txt=xlsread('//mnt/user02_3/trust5.xlsx');
trust_txt=readtable('//mnt/user02_3/trust5.xlsx','ReadVariableNames',1)
trust_txt([134 138 162],:)=[];
PMAT_CR=trust_txt(:,4);
trust_subj=table2array(trust_txt(:,1));
cd('//mnt/user02_3/structural/T1ImgGMV/');
for i=1:length(trust_subj)
    sub_num=num2str(trust_subj(i));
    if sub_num(1)=='1'
    if sub_num(2)=='0'
        sub_name=['S',sub_num(3:4)];
    else
        sub_name=['S',sub_num(2:4)];
    end
    else
        sub_name=['S',sub_num(1:4)];
    end
    name{i,1}=sub_name;
    cormatrix=dir(['*',sub_name,'_*']);
    infor=load_nifti(cormatrix.name)
    view_nii(infor)
end
 
%  mkdir('/mnt/user02_3/trust_MCI_264')
 cd('//mnt/user02_3/trust_MCI_268');
 save('CPM_gloCWF.mat','PMAT_CR','rest_1_mats');
 
 
 %%create and combine CPM map
 clc;clear;
% trust_txt=xlsread('//mnt/user02_3/trust5.xlsx');
trust_txt=readtable('//mnt/user02_3/trust6.xlsx','ReadVariableNames',1);
trust_txt([134 138 162],:)=[];
PMAT_CR=trust_txt(:,4);
trust_subj=table2array(trust_txt(:,1));
ggidx1=[find(trust_subj==5057)];
ggidx2=[find(trust_subj==5058)];


for i=1:length(trust_subj)
    if i<=ggidx1
cd('//mnt/user02_3/resting_analysis/Results160/ROISignals_FunImgARglobalCWF//');
    else
cd('//mnt/user02_3/resting_analysis2/Results160/ROISignals_FunImgARglobalCWF//');
    end
    sub_num=num2str(trust_subj(i));
    if sub_num(1)=='1'
    if sub_num(2)=='0'
        sub_name=['S',sub_num(3:4)];
    else
        sub_name=['S',sub_num(2:4)];
    end
    else
        sub_name=['S',sub_num(1:4)];
    end
    name{i,1}=sub_name;
    cormatrix=['ROICorrelation_FisherZ_',sub_name,'.mat'];
    load(cormatrix);
    rest_1_mats(:,:,i)=ROICorrelation_FisherZ;
end

%  mkdir('/mnt/user02_3/trust_MCI_264')
 cd('//mnt/user02_3/trust_MCI_160');
 save('CPM_gloCWF2.mat','PMAT_CR','rest_1_mats');
