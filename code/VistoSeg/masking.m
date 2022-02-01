function masking(filename, I)
Segmentations = load(filename);
O = fieldnames(Segmentations);

maskL = Segmentations.Lipofuscin;
maskM = Segmentations.MAP2;
maskG = Segmentations.GFAP;
Segmentations_m = Segmentations;
	 
for i = 1:numel(O)
	eval(['Segmentations_m.',O{i},'(maskL) = 0;']);
    eval(['Segmentations_m.',O{i},'(maskM) = 0;']);
    eval(['Segmentations_m.',O{i},'(maskG) = 0;']);
	disp(['Completed Masking ',O{i}]) 
end

Segmentations_m.DAPI = Segmentations.DAPI;

save([filename(1:end-4),'_masked.mat'],'-struct','Segmentations_m')

IMG = [];
for i = 1:numel(O)
IMG = [IMG, Segmentations_m.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
end
imwrite(IMG,[filename(1:end-4),'.png']);
end

