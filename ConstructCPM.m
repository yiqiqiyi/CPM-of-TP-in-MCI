clc;clear;
%load the subject ID
trust_txt=readtable('//mnt/user02_3/trust6_CPM.xlsx','ReadVariableNames',1)
%extract ID
trust_subj=table2array(trust_txt(:,1));
%CD to the folder with RSFC data
cd('//mnt/user02_3/resting_analysis2/Results160/ROISignals_FunImgARglobalCWF/');
%Construct the RSFC matrix based on the order of Subject ID.
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
 
 mkdir('/mnt/user02_3/trust_MCI_160')
 cd('//mnt/user02_3/trust_MCI_160');
 save('CPM_gloCWF2.mat','rest_1_mats');