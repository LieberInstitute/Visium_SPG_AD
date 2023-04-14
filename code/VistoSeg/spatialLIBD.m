function spatialLIBD(filename,filename1,slide,array,brain)

%filename = '/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Segmentations/VIFAD2_V10A27-106_B1_segmentation.mat';

n = [slide,'_',array,'_',brain];

load(filename)
[y,x] = size(pTau);
L = 600/x;
H = 2000/x;
lx = floor(L*x); ly = floor(L*y);
hx = ceil(H*x); hy = ceil(H*y);

AbetaL = imresize(Abeta,[ly lx]);  
AbetaH = imresize(Abeta,[hy hx]);
imwrite(AbetaL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_Abeta_seg_lowres.png'])
imwrite(AbetaH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_Abeta_seg_hires.png'])

pTauL = imresize(pTau,[ly lx]);  
pTauH = imresize(pTau,[hy hx]); 
imwrite(pTauL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_pTau_seg_lowres.png'])
imwrite(pTauH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_pTau_seg_hires.png'])

dapiL = imresize(DAPI,[ly lx]);  
dapiH = imresize(DAPI,[hy hx]); 
imwrite(dapiL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_DAPI_seg_lowres.png'])
imwrite(dapiH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_DAPI_seg_hires.png'])

mergeL = cat(3,pTauL,pTauL+AbetaL,dapiL);
mergeH = cat(3,pTauH,pTauH+AbetaH,dapiH);
imwrite(mergeL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_merge_seg_lowres.png'])
imwrite(mergeH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_merge_seg_hires.png'])

clear AbetaH AbetaL Abeta pTauH pTauL pTau DAPI dapiH dapiL
%filename1 = '/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1.mat';
load(filename1)

AbetaL = imresize(mat2gray(Abeta),[ly lx]);  
AbetaH = imresize(mat2gray(Abeta),[hy hx]); 
imwrite(AbetaL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_Abeta_lowres.png'])
imwrite(AbetaH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_Abeta_hires.png'])

pTauL = imresize(mat2gray(pTau),[ly lx]);  
pTauH = imresize(mat2gray(pTau),[hy hx]); 
imwrite(pTauL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_pTau_lowres.png'])
imwrite(pTauH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_pTau_hires.png'])

dapiL = imresize(mat2gray(DAPI),[ly lx]);  
dapiH = imresize(mat2gray(DAPI),[hy hx]); 
imwrite(dapiL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_DAPI_lowres.png'])
imwrite(dapiH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_DAPI_hires.png'])

mergeL = cat(3,pTauL,AbetaL,pTauL+dapiL);
mergeH = cat(3,pTauH,AbetaH,pTauH+dapiH);
imwrite(mergeL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_merge_lowres.png'])
imwrite(mergeH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_SPG_AD/processed-data/Images/spatialLIBD_images/',n,'_merge_hires.png'])


