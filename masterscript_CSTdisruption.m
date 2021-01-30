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
% Load FM scores
baselineFM=load(strcat(studydir, strokedir, 'baselineFM.mat'), 'basline')
baselineFM=baselineFM.basline;
finalFM= load(strcat(studydir, strokedir, 'finalFM.mat'), 'final')
finalFM=finalFM.final;
changeFM= load(strcat(studydir, strokedir, 'changeFM.mat'), 'change')
changeFM=changeFM.change;

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

%  atlasnames
names = readtable(strcat(studydir, 'processing/BrainStemAtlas/right_names.txt'));
atlas_names = table2array(names);

% Load flipped, binarized lesionMasks (in 1mm)
for i=1:23
    subnames{i}=strcat('SUB', num2str(i))
    lesion=read_avw(strcat(studydir,strokedir,'lesionMasks/right_flipped/SUB',num2str(i),'_lesion1mm_right.nii.gz')); %91x109x91
    lesionMask{i}=lesion;
end

% Load probabilistic brainstem atlases
for i=1:length(atlas_names)
    tract=read_avw(strcat(studydir,'processing/BrainStemAtlas/',cell2mat(atlas_names(i)),'_Atlas_swap.nii.gz')); %91x109x91
    tractfile_prob{i}=tract;
end

alltracts=loadTractNemo_rightonly_noalteration();
CST=alltracts(1:2);
CST_L=reshape(cell2mat(CST(1)), [1 902629]);
CST_R=reshape(cell2mat(CST(2)), [1 902629]);

% calculate "lesion load"
clear lesionload
for i=1:23
    lesion=lesionMask{i};
    for j=1
        atlas=tractfile_prob{j};
        lesionload{i,j}=calculate_lesion_load(lesion,atlas);
    end
end

ll=cell2mat(lesionload)

%proportional recovery
clear recovered
for i=1:23
    prop(i)=(100-baselineFM(i))*0.7+baselineFM(i)
    recovered(i)=finalFM(i)./prop(i)
end

bar(ll)
hold on
bar(recovered/1000)
plot(recovered,ll, '*')
