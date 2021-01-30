% July 17th 2020
% PCA of NeMo volumes, see if the PCs are similar?

disconnectivitydir='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/disconnectivity/NeMo2_outputs/july23_voxelwise_flippedBeforeNemo/'
% Load FM scores
baselineFM=load(strcat(studydir, strokedir, 'baselineFM.mat'), 'basline')
baselineFM=baselineFM.basline;
finalFM= load(strcat(studydir, strokedir, 'finalFM.mat'), 'final')
finalFM=finalFM.final;
changeFM= load(strcat(studydir, strokedir, 'changeFM.mat'), 'change')
changeFM=changeFM.change;
clear voxelNemo

cerebellumlabel=read_avw('/usr/share/fsl/5.0/data/atlases/MNI/MNI-maxprob-thr50-2mm.nii.gz');
cerebellumlabel(cerebellumlabel~=2)=1;
cerebellumlabel(cerebellumlabel==2)=0;

cortexlabel=read_avw('/usr/share/fsl/5.0/data/atlases/MNI/MNI-maxprob-thr0-2mm.nii.gz');

idx1=find(cortexlabel==3);% 4 lobes
idx2=find(cortexlabel==5);
idx3=find(cortexlabel==8);
idx4=find(cortexlabel==6);

cortex=[idx1;idx2;idx3;idx4];
notcortex=setdiff(1:902629,cortex);
cortexlabel(notcortex)=1;
cortexlabel(cortex)=0;
% Load voxelwise nemo.

for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
    nemo=read_avw(strcat(disconnectivitydir, 'SUB',num2str(i), '_lesion1mm_right_nemo_output_chacovol_res2mm_mean.nii.gz')); %91x109x91
    nemo_cortex=nemo.*~cortexlabel; %cortexlabel=0 where there is ctx, 1 elsewhere
   % nemo_cerebellum=nemo.*~cerebellumlabel; %cerebellumlabel=0 where there is cerebellum, 1 elsewhere
   % nemo=nemo_cortex+nemo_cerebellum;
    nemo=reshape(nemo_cortex,[1 902629]);
    voxelNemo{i}=nemo;
end

voxelN=cell2mat(voxelNemo');
[coeff,score,latent,tsquared,explaind]=pca(voxelN);

comp1=coeff(:,1);
comp2=coeff(:,2);
comp3=coeff(:,3);
comp4=coeff(:,4);
comp5=coeff(:,5);

rh_tracts_a=[4 6 8 10 12 13 15 17 19 21 23 25];
 nemodir ='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/BrainStemAtlas/NeMo/*'
a=dir(nemodir);
clear corel
for i=1:12
     r=rh_tracts_a(i)
     tract=read_avw(strcat(nemodir, a(r).name));
     disp(a(r).name)
    corel(i)=corr(comp1,reshape(tract,[902629 1]), 'Type', 'Pearson')
end
clear corl
for i=1
    r=rh_tracts_a(i)
     tract=read_avw(strcat(nemodir, a(r).name));
     disp(a(r).name);
    for j=1:100
        c=comp1(randperm(length(comp1)));
        corl{i,j}=corr(c,reshape(tract,[902629 1]), 'Type', 'Spearman');
    end
end
histogram(cell2mat(corl),100)
[rho,~]=corr(comp1,reshape(tract,[902629 1]), 'Type', 'Spearman')

%correlation between PC1 and each tract's ChaCo voxelwise scores.
bar(corel)
xticks(1:23)
xticklabels(atlas_names)


comp1_reshape=reshape(comp1(1:902629), [91 109 91]);
comp2_reshape=reshape(comp2(1:902629),[91 109 91]);
comp3_reshape=reshape(comp3(1:902629),[91 109 91]);
comp4_reshape=reshape(comp4(1:902629),[91 109 91]);
comp5_reshape=reshape(comp5(1:902629),[91 109 91]);

save_avw(comp1_reshape, strcat(disconnectivitydir,'component1.nii'),'f', [2 2 2 2]);
save_avw(comp2_reshape, strcat(disconnectivitydir,'component2.nii'),'f', [2 2 2 2]);
save_avw(comp3_reshape, strcat(disconnectivitydir,'component3.nii'),'f', [2 2 2 2]);
save_avw(comp4_reshape, strcat(disconnectivitydir,'component4.nii'),'f', [2 2 2 2]);
save_avw(comp5_reshape, strcat(disconnectivitydir,'component4.nii'),'f', [2 2 2 2]);

score1=score(:,1);
score2=score(:,2);
score3=score(:,3);
score4=score(:,4);
bar(vol)

bar(score1)
bar(score2)
bar(score3)
bar(score4);

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

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);
volnorm=vol/max(vol);
volnorm=volnorm*200;
%scatter(score1,score2, 'o', 'MarkerFaceColor', 'k')

for c=1:23
   scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   %scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end

xlabel('1st PC')
ylabel('2nd PC')
%text(score1(1:23), score2(1:23), subnames(1:23), 'FontSize', 10)
set(gca, 'FontSize', 13)

find(comp1>0.3)
% Explained variance
bar(explaind)
score1(15)=[];
baselineFM(15)=[];
%Correlation with baseline impairment
plot(score1,baselineFM, '*r')
%text(score1(1:22), baselineFM(1:22), subnames(1:22), 'FontSize', 8)
hold on;
[rho,p]=corrcoef(score1,baselineFM)
b=polyfit(score1, baselineFM,1);
a=polyval(b, score1);
plot(score1,a)


%Correlation with baseline impairment
plot(score2,baselineFM, '*r')
hold on;
[rho,p]=corrcoef(score2,baselineFM)
b=polyfit(score2, baselineFM,1);
a=polyval(b, score2);
plot(score2,a)

%Correlation with baseline impairment
plot(score3,baselineFM, '*r')
hold on;
[rho,p]=corrcoef(score3,baselineFM)
b=polyfit(score3, baselineFM,1);
a=polyval(b, score3);
plot(score3,a)


% input x rotated to new basis PCs
score1=score(:,1)
score2=score(:,2)

% Correlation with lesion volume.
plot3(score1,score2,vol,'*r')
grid on
xlabel('component1')
ylabel('component2')
text(score1(1:23), score2(1:23), vol(1:23), subnames(1:23), 'FontSize', 8)


plot(score1,vol, '*r')
text(score1(1:23), vol(1:23), subnames(1:23), 'FontSize', 8)
[rho,p]=corr(score1,vol)
