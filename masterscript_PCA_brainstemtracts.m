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

addpath([getenv('FSLDIR') '/etc/matlab']);
setenv( 'FSLDIR', '/usr/share/fsl/5.0');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldirmpath;

addpath([getenv('conn') '/etc/matlab']);
addpath([getenv('spm12') '/etc/matlab']);
studydir='/home/emo4002/colossus_shared3/pons_sfmodelling/'
resultsdir='results/ICC/'
strokedir='strokepts/'
controldir='control_subs/control_processed/'
conndir='preprocessing_conn/'

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

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))


%% Lesion overlap PCA
clear tractfile

% Load binarized brainstem atlases (in 1mm)
for i=1:length(atlas_names)
    tract=read_avw(strcat(studydir,'processing/BrainStemAtlas/',cell2mat(atlas_names(i)),'_Atlas_swap_thr0_1.nii.gz')); %91x109x91
    tractfile{i}=tract;
end

% Load probabilistic brainstem atlases
for i=1:length(atlas_names)
    tract=read_avw(strcat(studydir,'processing/BrainStemAtlas/',cell2mat(atlas_names(i)),'_Atlas_swap.nii.gz')); %91x109x91
    tractfile_prob{i}=tract;
end

% Load flipped, binarized lesionMasks (in 1mm)
for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
    lesion=read_avw(strcat(studydir,strokedir,'lesionMasks/right_flipped/SUB',num2str(i),'_lesion1mm_right.nii.gz')); %91x109x91
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

% calculate "lesion load"
clear lesionload
for i=1:23
    lesion=lesionMask{i};
    for j=1:12
        atlas=tractfile_prob{j};
        lesionload{i,j}=calculate_lesion_load(lesion,atlas);
    end
end


% simultaneous regression with PCA to account for lesion volume.
tune_lambda=acPCAtuneLambda(dc,lesionvol(:,1),0:0.1:10,11, 'linear',3,0)
lambda=tune_lambda.best_lambda;

outp = acPCA(dc,lesionvol(:,1),lambda,11, 'linear',3,1)

scores=outp.Xv; % projected data.
coeffs=outp.v; % principal components
figure(2)
bar(coeffs(:,1))
score1=scores(:,1)

figure(1)
plot(score1,baselineFM, '*r')
hold on;
[rho,p]=corr(score1,baselineFM)
b=polyfit(score1, baselineFM,1);
a=polyval(b, score1);
plot(score1,a)
title('AC-PC1 dice coefficient overlap all tracts vs. baseline FM')
xlabel('Dice overlap')
ylabel('Baseline FM')

figure(2)
plot(cst_resid,baselineFM, '*r')
hold on;
[rho,p]=corr(cst_resid,baselineFM)
b=polyfit(cst_resid, baselineFM,1);
a=polyval(b, cst_resid);
plot(cst_resid,a)
title('Residualized dice coefficient overlap with CST vs. baseline FM')
xlabel('Dice overlap')
ylabel('Baseline FM')

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
        Load_overlap{i,j}=calculate_lesion_load(subject_nemo',tract');
    end
end

sp=cell2mat(Load_overlap);
sp1=sp(1:14,:);
sp2=sp(16:23,:);
sp=[sp1;sp2];
% simultaneous regression with PCA to account for lesion volume.
tune_lambda=acPCAtuneLambda(sp,lesionvol(:,1),0:0.1:10,11, 'linear',3,0)
lambda=tune_lambda.best_lambda;

outp = acPCA(sp,lesionvol(:,1),lambda,11, 'linear',3,1)

score=outp.Xv; % projected data.
coeff=outp.v; % principal components

%Plot explained variance
bar(outp.varX_perc)

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

for c=1:22
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

clear tract_unaltered
tract_unaltered = loadTractNemo_rightonly_noalteration();
clear Load_overlap
for i=1:23
    subject_nemo=read_avw(strcat(disconnectivitydir, 'SUB',num2str(i), '_lesion1mm_nemo_output_chacovol_res2mm_mean_right.nii.gz'));
    subject_nemo=reshape(subject_nemo, [1 902629]);
    for j=1:12 % loop over tracts
        tract = tract_unaltered{j};
        tract=reshape(tract, [1 902629]);
        Load_overlap{i,j}=calculate_lesion_load(subject_nemo',tract');
    end
end

sp=cell2mat(Load_overlap);
sp1=sp(1:14,:);
sp2=sp(16:23,:);
sp=[sp1;sp2];

% simultaneous regression with PCA to account for lesion volume.
tune_lambda=acPCAtuneLambda(sp,lesionvol(:,1),0:0.1:10,11, 'linear',3,0)
lambda=tune_lambda.best_lambda;

outp = acPCA(sp,lesionvol(:,1),lambda,11, 'linear',3,1)

score=outp.Xv; % projected data.
coeff=outp.v; % principal components

%Plot explained variance
bar(outp.varX_perc)
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
bar(score2)
xticks(1:23)

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);

for c=1:22
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
