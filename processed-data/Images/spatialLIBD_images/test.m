filename = '/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Segmentations/VIFAD2_V10A27-106_B1_segmentation.mat';
filename1 = '/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1.mat';

load(filename)
[~,n,~] = fileparts(filename);
[y,x] = size(pTau);
L = 600/x;
H = 2000/x;
lx = floor(L*x); ly = floor(L*y);
hx = ceil(H*x); hy = ceil(H*y);

AbetaL = imresize(Abeta,[ly lx]);  
AbetaH = imresize(Abeta,[hy hx]); 
imwrite(AbetaL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'Abeta_seg_lowres.png'])
imwrite(AbetaH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'Abeta_seg_hires.png'])

pTauL = imresize(pTau,[ly lx]);  
pTauH = imresize(pTau,[hy hx]); 
imwrite(pTauL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'pTau_seg_lowres.png'])
imwrite(pTauH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'pTau_seg_hires.png'])

load(filename1)

AbetaL = imresize(Abeta,[ly lx]);  
AbetaH = imresize(Abeta,[hy hx]); 
imwrite(AbetaL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'Abeta_lowres.png'])
imwrite(AbetaH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'Abeta_hires.png'])

pTauL = imresize(pTau,[ly lx]);  
pTauH = imresize(pTau,[hy hx]); 
imwrite(pTauL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'pTau_lowres.png'])
imwrite(pTauH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'pTau_hires.png'])

dapiL = imresize(DAPI,[ly lx]);  
dapiH = imresize(DAPI,[hy hx]); 

mergeL = cat(3,pTauL,pTauL+AbetaL,zeros(ly,lx));
mergeH = cat(3,pTauH,pTauH+AbetaH,dapiH);
imwrite(AbetaL, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'Abeta_lowres.png'])
imwrite(AbetaH, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'Abeta_hires.png'])

pTauL = cat(3,pTau_hires,pTau_hires,zeros(ly,lx));
pTauH = cat(3,pTau_hires,pTau_hires,zeros(hy,hx));
imwrite(pTau_lowres, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'pTau_lowres.png'])
imwrite(pTau_hires, ['/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/spatialLIBD_images/',n(1:end-12),'pTau_hires.png'])

