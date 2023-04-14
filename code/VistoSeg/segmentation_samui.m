function segmentation_samui(fname)

load(fname);
img.DAPI = DAPI;
img.Abeta = Abeta;
img.pTau = pTau;
img.GFAP = GFAP;
img.MAP2 = MAP2;
img.Lipofuscin = Lipofuscin;

O = fieldnames(img);
%N = 4; %number of capture areas
disp(['The segmented image has ',num2str(numel(O)),' channels'])

imwrite(img.(O{1}),[fname(1:end-4),'.tif'])

for i = 2:numel(O)
imwrite(img.(O{i}),[fname(1:end-4),'.tif'],'writemode', 'append')
end
