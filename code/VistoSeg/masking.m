function masking(filename)
seg = load(filename);
O = fieldnames(seg);

Lip = contains(O,'Lipofuscin');  
mask = O{Lip};		 
eval(['mask = Segmentations.',mask,';']);
Segmentations_m = Segmentations;
	 
for i = 1:numel(O)
	 eval(['Segmentations_m.',O{i},'(mask) = 0;']);
	 disp(['Completed Masking ',O{i}]) 
end

save([filename(1:end-4),'_segmentation_masked.mat'],'Segmentations_m')
end

