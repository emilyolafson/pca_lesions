% Lesion location correlation with F-M outcome.


% atlasnames
names = readtable(strcat(studydir, 'processing/BrainStemAtlas/names.txt'));
atlas_names = table2array(names);

% Load binarized brainstem atlases (in 1mm)
for i=1:length(atlas_names)
    tract=read_avw(strcat(studydir,'processing/BrainStemAtlas/',cell2mat(atlas_names(i)),'_Atlas_swap_thr0_1.nii.gz')); %91x109x91
    tractfile{i}=tract;
end

% Load binarized lesionMasks (in 1mm)
for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
    lesion=read_avw(strcat(studydir,strokedir,'stroke_pts/lesionMasks/SUB',num2str(i),'_lesion1mm.nii.gz')); %91x109x91
    lesionMask{i}=lesion;
end

% Load FM scores
baselineFM=load(strcat(studydir, strokedir, 'baselineFM.mat'), 'basline')
baselineFM=baselineFM.basline;
finalFM= load(strcat(studydir, strokedir, 'finalFM.mat'), 'final')
finalFM=finalFM.final;
changeFM= load(strcat(studydir, strokedir, 'changeFM.mat'), 'change')
changeFM=changeFM.change;

% calculate dice coeff
for i=1:23
    lesion=lesionMask{i};
    for j=1:23
        atlas=tractfile{j};
        dice{i,j} = calculate_dice_coeff(lesion, atlas);
    end
end
p=1;
k=1;
clear dice_lr
for i=1:12
    if i==11
        dice_lr{p}=cell2mat(dice(:,k))'
        p=p+1
        k=k+1;
        continue
    end
    L_dice=dice(:,k)
    R_dice=dice(:,k+1)
    k=k+2;
    df = readmatrix(strcat(studydir,strokedir, 'demog_strokepts2.csv'));
    lr=df(:,2);
    left=find(lr==1);
    right=find(lr==0);
    clear a
    a(right)=R_dice(right);
    a(left)=L_dice(left);
    dice_lr{p}=cell2mat(a);
    p=p+1;
end

dice_lr=cell2mat(dice_lr')'
input_pca=[vol,dice_lr];
inp=normalize(input_pca);
idx=find(isnan(inp))
inp(idx)=0;

plot(cst,baselineFM, '*r');
[rho, p]=corr(cst',baselineFM)
hold on;
b=polyfit(cst', baselineFM,1);
a=polyval(b, cst');
plot(cst',a)

dc = cell2mat(dice);

fig1=figure(1)
set(fig1, 'Position', [ 300 300 1500 400])
bar(cell2mat(dice(1:23,:)'))
xticks(1:23)
xticklabels(atlas_names)
set(gca, 'FontSize', 10)
labels=subnames(1:23)';
legend(labels, 'Location', 'eastoutside');

input_pca=[vol, dc];
inpt=normalize(input_pca)
[coeff,score,latent,tsquared,explaind]=pca(dc);

% Explained variance across components
bar(explaind)

comp1=coeff(:,1)
comp2=coeff(:,2)
comp3=coeff(:,3)
comp4=coeff(:,4);
comp5=coeff(:,5);

bar(comp3)
bar(comp4)
atlas_names={'lesion vol', 'CST', 'FPT', 'ICPMC','ICPVC', 'LL', 'MCP','ML', 'POPT', 'SCPCR', 'SCPCT', 'SCPSC', 'STT'}
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

% principal component 3 - weights across variables
fig3=figure(3)
set(fig3, 'Position', [ 300 300 1500 400])
bar(comp3)
xticks(1:13)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 4 - weights across variables
fig4=figure(4)
set(fig4, 'Position', [ 300 300 1500 400])
bar(comp4)
xticks(1:13)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% principal component 5 - weights across variables
fig5=figure(5)
set(fig5, 'Position', [ 300 300 1500 400])
bar(comp5)
xticks(1:23)
xticklabels(atlas_names)
set(gca,'FontSize', 10)

% input x rotated to new basis PCs
score1=score(:,1)
score2=score(:,2)
plot3(score1,score2,vol,'*r')
grid on
xlabel('component1')
ylabel('component2')
text(score1(1:23), score2(1:23), vol(1:23), subnames(1:23), 'FontSize', 8)



score3=score(:,3)
plot(score3,baselineFM, '*r')
text(score3(1:23), baselineFM(1:23), subnames(1:23), 'FontSize', 8)

plot3(score1_a,score2_a,score3,'*r')

score4=score(:,4)

bar(score3)
bar(score4)

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
  %  scatter(score1(c),score2(c),volnorm(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   scatter(score1(c),score2(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
    hold on;
end
xlabel('1st PC')
ylabel('2nd PC')
%text(score1(1:23), score2(1:23), subnames(1:23), 'FontSize', 10)
set(gca, 'FontSize', 13)


plot(score2, vol, '*r')
hold on;
[rho,p]=corr(score2,vol)
b=polyfit(score2, vol,1);
a=polyval(b, score2);
plot(score2,a)


plot(score2,baselineFM, '*r')
hold on;
[rho,p]=corr(score2,baselineFM)
b=polyfit(score2, baselineFM,1);
a=polyval(b, score2);
plot(score2,a)

plot(score2,finalFM, '*r')
hold on;
[rho,p]=corr(score2,finalFM)
b=polyfit(score2, finalFM,1);
a=polyval(b, score2);
plot(score2,a)

plot(score2,changeFM, '*r')
hold on;
[rho,p]=corr(score2,changeFM)
b=polyfit(score2, changeFM,1);
a=polyval(b, score2);
plot(score2,a)

plot(vol,changeFM, '*r')
hold on;
[rho,p]=corr(vol,changeFM)
b=polyfit(vol, changeFM,1);
a=polyval(b, vol);
plot(vol,a)

for i=1:23
    c=base_z{i}% only cx + cb
    suma{i}=std(c)
end

basecorr=cell2mat(corr_degree_disconnect(:,1))
plot(score2,basecorr, '*r')
%[rho,p]=corr(score2, basecorr, 'Type', 'Spearman')
text(score2(1:23), basecorr(1:23), subnames(1:23), 'FontSize', 8)


plot(score2, cell2mat(suma),'*r')
[rho,p]=corr(cell2mat(suma)', score2)
bar(cell2mat(suma))
dice(2,:)

%subject scores PC2
bar(score2);
xticks(1:23)

bar(sumz');
xticks(1:23)

bstem_rois=[265,133,267,129,266,131,130,132,268,251,103,104,256,126,262];
cerebellum_rois=[100 102 114 116 113 118 109 107 106 105 117 115 110 243 242 247 246 241 252 245 238 201 254 253 237 239];
setdiff(bstem_rois,cerebellum_rois)
bstem_cereb=[bstem_rois, cerebellum_rois];
cortex = setdiff([1:268], bstem_cereb);

for i=1:23
    regional_disconnect = load(strcat(studydir, 'processing/disconnectivity/NeMo2_outputs/july13_shen/', 'SUB', num2str(i), '_lesion_1mmMNI_shen268_mean_chacovol.csv'));
   % regs=regional_disconnect(setdiff(1:268,bstem_rois));
    sumchaco_notbrainstem{i} = sum(regs)
    sum_cereb=sum(regional_disconnect(cerebellum_rois))
    sum_cortex=sum(regional_disconnect(cortex))
    wtd_sum{i}=sum_cereb*(1/7)+sum_cortex*(6/7);
end

sumz=cell2mat(wtd_sum)


plot(sumz',score2, '*r')
hold on;
[rho,p]=corr(sumz',score2)
b=polyfit(sumz', score2,1);
a=polyval(b, sumz');
plot(sumz',a)
text(sumz(1:23), score2(1:23), subnames(1:23), 'FontSize', 8)


csvwrite(strcat(studydir,strokedir,'dice_brainstematlases.csv'), dice)



%% ICA

obj=rica(dc,4)
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

plot(score1,score2, '*r')


for c=1:23
 %  scatter(score2(c),score4(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
  plot3(score1(c),score2(c),baselineFM(c),[colorstring(c), 'o'], 'MarkerFaceColor', colorstring(c))
   xlabel('1st PC')
  ylabel('2nd PC')
   grid on
    hold on;
end

%% Multimodality?
plot(dc(:,1))

