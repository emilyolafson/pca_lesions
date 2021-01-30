% July 17th 2020
% PCA of NeMo volumes, see if the PCs are similar?

disconnectivitydir='processing/disconnectivity/NeMo2_outputs/july17_voxelwise/'
% Load FM scores
baselineFM=load(strcat(studydir, strokedir, 'baselineFM.mat'), 'basline')
baselineFM=baselineFM.basline;
finalFM= load(strcat(studydir, strokedir, 'finalFM.mat'), 'final')
finalFM=finalFM.final;
changeFM= load(strcat(studydir, strokedir, 'changeFM.mat'), 'change')
changeFM=changeFM.change;

% Load voxelwise nemo.
for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
   % nemo=read_avw(strcat(studydir,disconnectivitydir,'/SUB',num2str(i),'_lesion1mm_nemo_output_chacovol_res2mm_smooth6mm_mean.nii.gz')); %91x109x91
   % nemo=reshape(nemo,[1 902629]);
   % voxelNemo{i}=nemo;
end

voxelN=cell2mat(voxelNemo');
lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
vol=lesionvol(:,2);
input_pca = [vol,voxelN];
input_pca=normalize(input_pca);
idx=find(isnan(input_pca));
input_pca(idx)=0;
[coeff,score,latent,tsquared,explaind]=pca(input_pca);

comp1=coeff(:,1);
comp2=coeff(:,2);
comp3=coeff(:,3);
comp4=coeff(:,4);

comp1_reshape=reshape(comp1(1:902629), [91 109 91]);
comp2_reshape=reshape(comp2(1:902629),[91 109 91]);
comp3_reshape=reshape(comp3(1:902629),[91 109 91]);
comp_reshape=reshape(comp4(1:902629),[91 109 91]);

save_nii(make_nii(comp1_reshape, [2 2 2]), strcat(disconnectivitydir,'component1.nii'));
save_nii(make_nii(comp2_reshape, [2 2 2]), strcat(disconnectivitydir,'component2.nii')); 
save_nii(make_nii(comp2_reshape, [2 2 2]), strcat(disconnectivitydir,'component3.nii'));
save_nii(make_nii(comp2_reshape, [2 2 2]), strcat(disconnectivitydir,'component4.nii'));

score1=score(:,1);
score2=score(:,2);
score3=score(:,3);
score4=score(:,4);

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
   scatter(score1(c),score2(c),volnorm(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   %scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end
xlabel('1st PC')
ylabel('2nd PC')
%text(score1(1:23), score2(1:23), subnames(1:23), 'FontSize', 10)
set(gca, 'FontSize', 13)


% Explained variance
bar(explaind)

%Correlation with baseline impairment
plot(score1,baselineFM, '*r')
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
