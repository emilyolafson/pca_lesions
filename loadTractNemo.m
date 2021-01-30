function [tractz] = loadTractNemo()
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
   
    for i = 1:23
        z=i+2;
        tract=read_avw(strcat(nemodir, a(z).name));
        if (i==5 || i==6 || i==7 || i==8 || i==9 || i==11 || i > 14) % cerebellum labels
            disp([num2str(i), 'one'])
            tract=tract.*cortexlabel; 
        end
        
        if (i < 5 || i==10 || i==12 || i==13 || i==14 || i==15) % cortical labels
            disp([num2str(i), 'two'])
            tract=tract.*cerebellumlabel;
        end
        tractz{i}=tract;
        save_avw(tractz{i},strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/processing/BrainStemAtlas/altered/', a(z).name), 'f', [2 2 2 2])
    end
end

