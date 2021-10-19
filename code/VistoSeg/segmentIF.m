function segmentIF(filename,I)

img = load(filename);
O = fieldnames(img);

for i = 1:numel(O)
    Re.(O{i}) = rescale(img.(O{i}));
    T.(O{i}) = graythresh(Re.(O{i}));  
    B.(O{i}) = imbinarize(Re.(O{i}),T.(O{i})); 
end

M = B;

figure('Visible','off')
for i = 1:numel(O)
subplot(2,3,i)
imhist(Re.(O{i})(:));
ylim([0 2.5*10^7])
xline(T.(O{i}),'r',{num2str(T.(O{i}))});
title(O{i})
M.(O{i})(B.Lipofuscin|B.MAP2|B.GFAP) = 0;
end
saveas(gcf,[filename(1:end-4),'_thresh.png'])
close all

IMG = [];
IMG1 = [];
IMG2 = [];
for i = 1:numel(O)
IMG = [IMG,Re.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
IMG1 = [IMG1,B.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
IMG2 = [IMG2,M.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
end

BW = bwareaopen(M.Abeta, 200);
M.Abeta = BW;
M.DAPI = B.DAPI;
IMG3 = [];
for i = 1:numel(O)
IMG3 = [IMG3,M.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
end

imwrite([IMG;ones(20,size(IMG,2));IMG1;ones(20,size(IMG,2));IMG2;ones(20,size(IMG,2));IMG3],[filename(1:end-4),'_test.png'])
imwrite(BW,[filename(1:end-4),'_Abeta.png'])
imwrite(M.pTau,[filename(1:end-4),'_pTau.png'])

save([filename(1:end-4),'_segmentation.mat'], '-struct', 'M')		 

end