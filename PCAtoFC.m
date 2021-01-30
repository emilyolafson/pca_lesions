% Bridge PCA to FC
save(strcat(studydir,resultsdir, 'base_degree_allsubjects.mat'), 'base_z')
save(strcat(studydir,resultsdir, 'base_vs_FU_change_degree_allsubjects.mat'), 'change_degree')

studydir='/home/emo4002/colossus_shared3/pons_sfmodelling/'
resultsdir='results/ICC/'
cortexlabel=read_avw('/usr/share/fsl/5.0/data/atlases/MNI/MNI-maxprob-thr0-2mm.nii.gz');

idx1=find(cortexlabel==3);% 4 lobes
idx2=find(cortexlabel==5);
idx3=find(cortexlabel==8);
idx4=find(cortexlabel==6);

cortex=[idx1;idx2;idx3;idx4];
notcortex=setdiff(1:902629,cortex);
cortexlabel(notcortex)=1;
cortexlabel(cortex)=0;

comp1_right=read_avw(strcat(disconnectivitydir,'component1_right.nii'));
%imagesc(squeeze(comp1_right(:,30,:)))
comp1_right=reshape(comp1_right,[902629 1]);
comp1_left=read_avw(strcat(disconnectivitydir,'component1_left.nii'));
%imagesc(squeeze(comp1_left(:,30,:)))
comp1_left=reshape(comp1_left,[902629 1]);

clear icc

df = readmatrix(strcat(studydir,strokedir, 'demog_strokepts2.csv'));
lr=df(:,2);
left=find(lr==1);
right=find(lr==0);

for i=1:23
    j=3
    if i==6
    %   j=4
    end
    if i==12
      %  j=3
    end
    if i==20
        %
       j=2
    end
    ICC_sub=read_avw(strcat(studydir,'results/ICC/patients/icc_makenii/SUB',num2str(i), '_S',num2str(j), '_ICC_swap.nii.gz'));
    %bold_sub=read_avw(strcat(studydir,resultsdir,'SUB', num2str(i),'_S1_meanbold_GSR.nii.gz'));
    if(ismember(i,right))
       a=reshape(ICC_sub, [902629 1]);
        icc{i}=a.*comp1_right;

        correlz{i}=corr(a(comp1_right>0),comp1_right(comp1_right>0), 'Type', 'Spearman')
    end
    if(ismember(i,left))
       a=reshape(ICC_sub, [902629 1]);
        icc{i}=a.*comp1_left;
        [rho,p]=corr(a(comp1_left>0),comp1_left(comp1_left>0), 'Type', 'Spearman')
        scatter(log(a(comp1_left>0)),log(comp1_left(comp1_left>0)), 'o', 'Filled','MarkerFaceAlpha', 0.0020)
       correlz{i}=corr(a(comp1_left>0),comp1_left(comp1_left>0), 'Type', 'Spearman')

    end
end

plot(score1, cell2mat(correlz)','*r')
[rho,p]=corr(cell2mat(correlz)',score1)

for i=1:23
    meanicc{i}=mean(icc{i});
end
meanz=cell2mat(meanicc);

df = readmatrix(strcat(studydir,strokedir, 'demog_strokepts2.csv'));
lr=df(:,2);
left=find(lr==1);
right=find(lr==0);
colorstring='';
for i=1:size(left,1)
    colorstring(left(i))='b';
end
for i=1:size(right,1)
    colorstring(right(i))='r';
end
clf
scatter(score1,meanz, '*r')
hold on;
[rho,p]=corr(score1,meanz', 'Type', 'Pearson')
b=polyfit(score1, meanz',1);
a=polyval(b, score1);
plot(score1,a)
title('S5')

for c=1:23
   scatter(score1(c),meanz(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   %scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end

text(score1(1:23), meanz(1:23), subnames(1:23), 'FontSize', 10)

plot(score1,vol, '*r')

bar(score1)
[rho,p]=corr(score1,meanz', 'Type', 'Pearson')


