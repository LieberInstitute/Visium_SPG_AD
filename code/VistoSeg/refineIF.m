function refineIF(filename,I,a,g,m,p)
img = load(filename);
O = fieldnames(img);

for i = 1:numel(O)
    Re.(O{i}) = rescale(img.(O{i}));
    T.(O{i}) = graythresh(Re.(O{i})); 
end
clear img

figure('Visible','off')
for i = 1:numel(O)
subplot(2,3,i)
imhist(Re.(O{i})(:));
ylim([0 2.5*10^7])
xline(T.(O{i}),'r',{num2str(T.(O{i}))});

if contains('Abeta',O{i}) && ~isempty(a)
    T.Abeta = a; 
    xline(T.(O{i}),'g',{num2str(T.(O{i}))});
elseif contains('GFAP',O{i}) && ~isempty(g)
    T.GFAP = g; 
    xline(T.(O{i}),'g',{num2str(T.(O{i}))});
elseif contains('MAP2',O{i}) && ~isempty(m)
    T.MAP2 = m; 
    xline(T.(O{i}),'g',{num2str(T.(O{i}))});
elseif contains('pTau',O{i}) && ~isempty(p)
    T.pTau = p; 
    xline(T.(O{i}),'g',{num2str(T.(O{i}))});
end
title(O{i})

B.(O{i}) = imbinarize(Re.(O{i}),T.(O{i})); 
end
saveas(gcf,[filename(1:end-4),'_thresh.png'])
close all

M = B;
for i = 1:numel(O)
M.(O{i})(B.Lipofuscin|B.MAP2|B.GFAP) = 0;
end


IMG = [];
IMG1 = [];
IMG2 = [];
for i = 1:numel(O)
IMG = [IMG,Re.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
IMG1 = [IMG1,B.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
IMG2 = [IMG2,M.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
end

BW = bwareaopen(M.Abeta, 1500);
%M.Abeta = BW;
M.Abeta = bwpropfilt(BW, 'Solidity', [0.5, 1]); 
M.DAPI = B.DAPI;
IMG3 = [];
for i = 1:numel(O)
IMG3 = [IMG3,M.(O{i})(I(1):I(2),I(3):I(4)),ones(1000,20)];
end

imwrite([IMG;ones(20,size(IMG,2));IMG1;ones(20,size(IMG,2));IMG2;ones(20,size(IMG,2));IMG3],[filename(1:end-4),'_test.png'])
imwrite(M.Abeta,[filename(1:end-4),'_Abeta.png'])
imwrite(M.pTau,[filename(1:end-4),'_pTau.png'])

save([filename(1:end-4),'_segmentation.mat'], '-struct', 'M')

end
