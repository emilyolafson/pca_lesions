function [tractz] = loadTractNemo_rightonly()
% Load voxelwise ChaCo scores for each tract & edit with appropriate mask.
    nemodir ='/home/emo4002/colossus_shared3/pons_sfmodelling/processing/BrainStemAtlas/NeMo/*'
    a=dir(nemodir);
    tractnames=a.name;
    

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
    
    rh_tracts=[2 4 6 8 10 11 13 15 17 19 21 23];
    rh_tracts_a=[4 6 8 10 12 13 15 17 19 21 23 25];
   
    for i = 1:size(rh_tracts,2)
        r=rh_tracts_a(i)
        tract=read_avw(strcat(nemodir, a(r).name));
        if (r==8 || r==10 || r==13 || r>17) % cerebellum labels
            disp([num2str(r), 'one'])
            tract=tract.*cortexlabel; 
        end
        
        if (r < 7 || r==12 || r==15 || r==17) % cortical labels
            disp([num2str(r), 'two'])
            tract=tract.*cerebellumlabel;
        end
        tractz{i}=tract;
        save_avw(tractz{i},strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/processing/BrainStemAtlas/', a(r).name), 'f', [2 2 2 2])
    end
end

