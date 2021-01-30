% Calculate the Spearman rank correlation between each subject's voxelwise
% ChaCo map and each tract's map.

addpath([getenv('FSLDIR') '/etc/matlab']);
setenv( 'FSLDIR', '/usr/share/fsl/5.0');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldirmpath;
studydir='/home/emo4002/colossus_shared3/pons_sfmodelling/'
resultsdir='results/degreeFC/'
strokedir='strokepts/'
controldir='control_subs/'
conndir='preprocessing_conn/'
figuresdir = 'results/figures/degree/';
disconnectivitydir = 'processing/disconnectivity/NeMo2_outputs/july17_voxelwise/';

baselineFM=load(strcat(studydir, strokedir, 'baselineFM.mat'), 'basline')
baselineFM=baselineFM.basline;


tractz = loadTractNemo();

for i=1:23
    subject_nemo=read_avw(strcat(studydir,disconnectivitydir, 'SUB',num2str(i), '_lesion1mm_nemo_output_chacovol_res2mm_mean.nii.gz'));
    subject_nemo=reshape(subject_nemo, [1 902629]);
    for j=1:23 % loop over tracts
        tract = tractz{j};
        tract=reshape(tract, [1 902629]);
        Spearman_overlap{i,j}=corr(subject_nemo',tract', 'Type', 'Spearman');
    end
end
Sp=cell2mat(Spearman_overlap);

plot(Sp(:,1))
xticks(1:23)
xticklabels(atlas_names)

[coeff,score,latent,tsquared,explaind]=pca(cell2mat(Spearman_overlap));

%Plot explained variance
bar(explaind)

comp1=coeff(:,1);
comp2=coeff(:,2);
comp3=coeff(:,3);
comp4=coeff(:,4);

% atlasnames
names = readtable(strcat(studydir, 'processing/BrainStemAtlas/names.txt'));
atlas_names = table2array(names);

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

%% PLOT 3D with labels & lines of best fit.
lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
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
% Load voxelwise nemo.
for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
end


fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);
volnorm=vol/max(vol);
volnorm=volnorm*200;
%scatter(score1,score2, 'o', 'MarkerFaceColor', 'k')
for c=1:23
   scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
  %plot3(score1(c),score2(c),baselineFM(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   xlabel('1st PC')
  ylabel('2nd PC')
   grid on
    hold on;
end
xyz=[score1(right),score2(right),baselineFM(right)];
r=mean(xyz,1);
xyz=bsxfun(@minus,xyz,r);
[~,~,V]=svd(xyz,0);
x_fit=@(z_fit) r(1)+(z_fit-r(3))/V(3,1)*V(1,1);
y_fit=@(z_fit) r(2)+(z_fit-r(3))/V(3,1)*V(2,1);
plot3(x_fit(baselineFM(right)),y_fit(baselineFM(right)),baselineFM(right),'r')
xyz=[score1(left),score2(left),baselineFM(left)];
r=mean(xyz,1);
xyz=bsxfun(@minus,xyz,r);
[~,~,V]=svd(xyz,0);
x_fit=@(z_fit) r(1)+(z_fit-r(3))/V(3,1)*V(1,1);
y_fit=@(z_fit) r(2)+(z_fit-r(3))/V(3,1)*V(2,1);
plot3(x_fit(baselineFM(left)),y_fit(baselineFM(left)),baselineFM(left),'b')
text(score1(1:23), score2(1:23), baselineFM(1:23), subnames(1:23), 'FontSize', 8)
grid on


pc1pc2=score1+score2
plot(pc1pc2,baselineFM, '*r')
[rho,p]=corrcoef(pc1pc2,baselineFM)
b=polyfit(pc1pc2, baselineFM,1);
a=polyval(b, pc1pc2);
hold on;
plot(pc1pc2,a)

%% ICA analysis; not assuming linearity
obj=rica(cell2mat(Spearman_overlap),4)
a=obj.TransformWeights;

z=corr(cell2mat(Spearman_overlap), coeff1)
bar(z)

coeff1=a(:,1)
coeff2=a(:,2)
coeff3=a(:,3)
coeff4=a(:,4)

fig1=figure(1)
set(fig1, 'Position', [ 300 300 1500 400])
bar(coeff1)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 2 - weights across variables
fig2=figure(2)
set(fig2, 'Position', [ 300 300 1500 400])
bar(coeff2)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 3 - weights across variables
fig3=figure(3)
set(fig3, 'Position', [ 300 300 1500 400])
bar(coeff3)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 2 - weights across variables
fig2=figure(2)
set(fig2, 'Position', [ 300 300 1500 400])
bar(coeff4)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)


scores=transform(obj,cell2mat(Spearman_overlap))

score1=scores(:,1)
score2=scores(:,2)
score3=scores(:,3)
score4=scores(:,4)

bar(score1)

plot(score3,score2, '*r')

for c=1:23
  scatter(score1(c),score3(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
  %plot3(score1(c),score2(c),baselineFM(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   xlabel('1st PC')
  ylabel('2nd PC')
   grid on
    hold on;
end

plot(score4,baselineFM, '*r')







names_noLR={'CST', 'FPT', 'ICPMC', 'ICPVC', 'LL', 'MCP', 'ML', 'POTPT', 'SCPCR', 'SCPCT', 'SCPSC', 'STT'};

% principal component 1 - weights across variables
fig1=figure(1)
set(fig1, 'Position', [ 300 300 1500 400])
bar(comp1)
xticks(1:12)
xticklabels(names_noLR)
set(gca,'FontSize', 10)

% principal component 2 - weights across variables
fig2=figure(2)
set(fig2, 'Position', [ 300 300 1500 400])
bar(comp2)
xticks(1:12)
xticklabels(names_noLR)
set(gca,'FontSize', 10)

lesionvol=load(strcat(studydir, 'strokepts/allpts_lesionvol.txt'))
fig3=figure(3)
set(fig3, 'Position', [0 0 700 700])
vol=lesionvol(:,2);
volnorm=vol/max(vol);
volnorm=volnorm*200;
%scatter(score1,score2, 'o', 'MarkerFaceColor', 'k')
for c=1:23
  % scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   plot3(score1(c),score2(c),baselineFM(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   xlabel('1st PC')
ylabel('2nd PC')
   grid on
    hold on;
end

