% July 24 2020

% flip lesions to all be on one side to increase power.
% Lesion PCA: Dice coefficient between tracts & lesion.
    % No parameters to change here.
% Tract PCA: 
    % Run this twice. Once calculating the Spearman correlation between
    % atlas tracts and lesion tracts flipped L -> R *after* NeMo2, one with
    % the lesion flipped L -> R *before* NeMo2.
    % Comparison between the two outputs: quantifyBeforeAfterChaCoflip.m
%% globals
%  atlasnames
names = readtable(strcat(studydir, 'processing/BrainStemAtlas/right_names.txt'));
atlas_names = table2array(names);

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

% Load FM scores
baselineFM=load(strcat(studydir, strokedir, 'baselineFM.mat'), 'basline')
baselineFM=baselineFM.basline;
finalFM= load(strcat(studydir, strokedir, 'finalFM.mat'), 'final')
finalFM=finalFM.final;
changeFM= load(strcat(studydir, strokedir, 'changeFM.mat'), 'change')
changeFM=changeFM.change;

%% Lesion overlap PCA
clear tractfile
% Load binarized brainstem atlases (in 1mm)
for i=1:length(atlas_names)
    tract=read_avw(strcat(studydir,'processing/BrainStemAtlas/',cell2mat(atlas_names(i)),'_Atlas_swap_thr0_1.nii.gz')); %91x109x91
    tractfile{i}=tract;
end

% Load flipped, binarized lesionMasks (in 1mm)
for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
    lesion=read_avw(strcat(studydir,strokedir,'stroke_pts/lesionMasks/right_flipped/SUB',num2str(i),'_lesion1mm_right.nii.gz')); %91x109x91
    lesionMask{i}=lesion;
end

% calculate dice coeff
clear dice
for i=1:23
    lesion=lesionMask{i};
    for j=1:12 % only calculate Dice with RIGHT atlas tracts + MCP
        atlas=tractfile{j};
        dice{i,j} = calculate_dice_coeff(lesion, atlas);
    end
end

dc=cell2mat(dice);
[coeff,score,latent,tsquared,explaind]=pca(dc);

% Explained variance across components
bar(explaind)

comp1=coeff(:,1)
comp2=coeff(:,2)
comp3=coeff(:,3)
comp4=coeff(:,4);
comp5=coeff(:,5);

% input x rotated to new basis PCs
score1=score(:,1)
score2=score(:,2)

bar(score1)
xticks(1:23)
set(gca,'FontSize', 10)

atlas_names={'CST', 'FPT', 'ICPMC','ICPVC', 'LL', 'MCP','ML', 'POPT', 'SCPCR', 'SCPCT', 'SCPSC', 'STT'}

% principal component 1 - weights across variables
fig1=figure(1)
set(fig1, 'Position', [ 300 300 1500 400])
bar(comp1)
xticks(1:13)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 2 - weights across variables
fig2=figure(2)
set(fig2, 'Position', [ 300 300 1500 400])
bar(comp2)
xticks(1:13)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);
volnorm=vol/max(vol);
volnorm=volnorm*200;
%scatter(score1,score2, 'o', 'MarkerFaceColor', 'k')
for c=1:23
  %  scatter(score1(c),score2(c),volnorm(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end
xlabel('1st PC')
ylabel('2nd PC')
%text(score1(1:23), score2(1:23), subnames(1:23), 'FontSize', 10)
set(gca, 'FontSize', 13)

plot(score1,vol,'*r')
[rho,p]=corr(score1,vol)

plot(score1,baselineFM, '*r')
hold on;
[rho,p]=corr(score1,baselineFM)
b=polyfit(score1, baselineFM,1);
a=polyval(b, score1);
plot(score1,a)

%% ChaCo scores overlap PCA
% flipped before NeMo 
disconnectivitydir='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/disconnectivity/NeMo2_outputs/july23_voxelwise_flippedBeforeNemo/'

tractz = loadTractNemo_rightonly();
clear Spearman_overlap
for i=1:23
    subject_nemo=read_avw(strcat(disconnectivitydir, 'SUB',num2str(i), '_lesion1mm_right_nemo_output_chacovol_res2mm_mean.nii.gz'));
    subject_nemo=reshape(subject_nemo, [1 902629]);
    for j=1:12 % loop over tracts
        tract = tractz{j};
        tract=reshape(tract, [1 902629]);
        Spearman_overlap{i,j}=corr(subject_nemo',tract', 'Type', 'Spearman');
    end
end

sp=cell2mat(Spearman_overlap);

[coeff,score,latent,tsquared,explaind]=pca(sp);

%Plot explained variance
bar(explaind)

comp1=coeff(:,1);
comp2=coeff(:,2);
comp3=coeff(:,3);
comp4=coeff(:,4);


% principal component 1 - weights across variables
fig1=figure(1)
set(fig1, 'Position', [ 300 300 1500 400])
bar(comp1)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 2 - weights across variables
fig2=figure(2)
set(fig2, 'Position', [ 300 300 1500 400])
bar(comp2)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 3 - weights across variables
fig3=figure(3)
set(fig3, 'Position', [ 300 300 1500 400])
bar(comp3)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

score1=score(:,1);
score2=score(:,2);
score3=score(:,3);

%plot across subjects
bar(score1)
xticks(1:23)

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);

for c=1:23
  %  scatter(score1(c),score2(c),volnorm(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end

plot(score1,baselineFM, '*r')
hold on;
[rho,p]=corr(score1,baselineFM)
b=polyfit(score1, baselineFM,1);
a=polyval(b, score1);
plot(score1,a)

clf
plot(score2,baselineFM, '*r')
hold on;
[rho,p]=corr(score2,baselineFM)
b=polyfit(score2, baselineFM,1);
a=polyval(b, score2);
plot(score2,a)


plot(score2,vol, '*r')
hold on;
[rho,p]=corr(score2,vol)
b=polyfit(score2, vol,1);
a=polyval(b, score2);
plot(score2,a)


% flipped after NeMo

disconnectivitydir='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/disconnectivity/NeMo2_outputs/july20_voxelwise_flippedAfterNeMo/'

clear tractz
tractz = loadTractNemo_rightonly();
clear Spearman_overlap
for i=1:23
    subject_nemo=read_avw(strcat(disconnectivitydir, 'SUB',num2str(i), '_lesion1mm_nemo_output_chacovol_res2mm_mean_right.nii.gz'));
    subject_nemo=reshape(subject_nemo, [1 902629]);
    for j=1:12 % loop over tracts
        tract = tractz{j};
        tract=reshape(tract, [1 902629]);
        Spearman_overlap{i,j}=corr(subject_nemo',tract', 'Type', 'Spearman');
    end
end

sp=cell2mat(Spearman_overlap);

[coeff,score,latent,tsquared,explaind]=pca(sp);

%Plot explained variance
bar(explaind)

comp1=coeff(:,1);
comp2=coeff(:,2);
comp3=coeff(:,3);
comp4=coeff(:,4);

sp
% principal component 1 - weights across variables
fig1=figure(1)
set(fig1, 'Position', [ 300 300 1500 400])
bar(comp1)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 2 - weights across variables
fig2=figure(2)
set(fig2, 'Position', [ 300 300 1500 400])
bar(comp2)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 3 - weights across variables
fig3=figure(3)
set(fig3, 'Position', [ 300 300 1500 400])
bar(comp3)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

score1=score(:,1);
score2=score(:,2);
score3=score(:,3);

%plot across subjects
bar(score2)
xticks(1:23)

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);

for c=1:23
   %plot3(score1(c),score2(c),vol(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end

plot(score1,baselineFM, '*r')
hold on;
[rho,p]=corr(score1,baselineFM)
b=polyfit(score1, baselineFM,1);
a=polyval(b, score1);
plot(score1,a)

clf
plot(score2,baselineFM, '*r')
hold on;
[rho,p]=corr(score2,baselineFM)
b=polyfit(score2, baselineFM,1);
a=polyval(b, score2);
plot(score2,a)

clf
plot(score2,vol, '*r')
hold on;
[rho,p]=corr(score2,vol)
b=polyfit(score2, vol,1);
a=polyval(b, score2);
plot(score2,a)
