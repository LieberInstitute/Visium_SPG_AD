function segmentIF(filename,I)

img = load(filename);
DAPI = rescale(img.DAPI);
Abeta = rescale(img.Abeta);
MAP2 = rescale(img.MAP2);
GFAP = rescale(img.GFAP);
LIP = rescale(img.Lipofuscin);
pTau = rescale(img.pTau);

figure('Visible','off')
[counts,x] = imhist(Abeta(:));
stem(x,counts,'g')
hold on
[counts,x] = imhist(MAP2(:));
stem(x,counts,'m')
[counts,x] = imhist(GFAP(:));
stem(x,counts,'r')
[counts,x] = imhist(LIP(:));
stem(x,counts,'c')
[counts,x] = imhist(DAPI(:));
stem(x,counts,'b')
[counts,x] = imhist(pTau(:));
stem(x,counts,'y')
ylim([0 2.5*10^7])

saveas(gcf,[filename(1:end-4),'_hists.png'])
close all

TDAPI = graythresh(DAPI);
TAbeta = graythresh(Abeta);
TMAP2 = graythresh(MAP2);
TGFAP = graythresh(GFAP);
TLIP = graythresh(LIP);
TpTau = graythresh(pTau);

figure('Visible','off')
subplot(2,3,1)
imhist(Abeta(:));
ylim([0 2.5*10^7])
xline(TAbeta,'r')
title(['Abeta ',num2str(TAbeta),''])

subplot(2,3,2)
imhist(MAP2(:));
ylim([0 2.5*10^7])
xline(TMAP2,'r')
title(['MAP2 ',num2str(TMAP2),''])

subplot(2,3,3)
imhist(GFAP(:));
ylim([0 2.5*10^7])
xline(TGFAP,'r')
title(['GFAP ',num2str(TGFAP),''])

subplot(2,3,4)
imhist(LIP(:));
ylim([0 2.5*10^7])
xline(TLIP,'r')
title(['Lip ',num2str(TLIP),''])

subplot(2,3,5)
imhist(DAPI(:));
ylim([0 2.5*10^7])
xline(TDAPI,'r')
title(['DAPI ',num2str(TDAPI)])

subplot(2,3,6)
imhist(pTau(:));
ylim([0 2.5*10^7])
xline(TpTau,'r')
title(['pTau ',num2str(TpTau),''])

saveas(gcf,[filename(1:end-4),'_thresh.png'])
close all

BDAPI = imbinarize(DAPI,TDAPI);
BLIP = imbinarize(LIP,TLIP);

BAbeta = imbinarize(Abeta,TAbeta);
BMAP2 = imbinarize(MAP2,TMAP2);
BGFAP = imbinarize(GFAP,TGFAP);
BpTau = imbinarize(pTau,TpTau);

O = fieldnames(img);
IMG = [Abeta(I(1):I(2),I(3):I(4)),ones(1000,20),MAP2(I(1):I(2),I(3):I(4)),ones(1000,20),mat2gray(GFAP(I(1):I(2),I(3):I(4))),ones(1000,20),LIP(I(1):I(2),I(3):I(4)),ones(1000,20),DAPI(I(1):I(2),I(3):I(4)),ones(1000,20),pTau(I(1):I(2),I(3):I(4))];
IMG1 = [BAbeta(I(1):I(2),I(3):I(4)),ones(1000,20),BMAP2(I(1):I(2),I(3):I(4)),ones(1000,20),BGFAP(I(1):I(2),I(3):I(4)),ones(1000,20),BLIP(I(1):I(2),I(3):I(4)),ones(1000,20),BDAPI(I(1):I(2),I(3):I(4)),ones(1000,20),BpTau(I(1):I(2),I(3):I(4))];

BAbeta(BLIP|BMAP2|BGFAP) = 0;
BpTau(BLIP|BMAP2|BGFAP) = 0;
BGFAP(BLIP|BMAP2|BGFAP) = 0;
BMAP2(BLIP|BMAP2|BGFAP) = 0;
BLIP(BLIP|BMAP2|BGFAP) = 0;

IMG2 = [BAbeta(I(1):I(2),I(3):I(4)),ones(1000,20),BMAP2(I(1):I(2),I(3):I(4)),ones(1000,20),BGFAP(I(1):I(2),I(3):I(4)),ones(1000,20),BLIP(I(1):I(2),I(3):I(4)),ones(1000,20),BDAPI(I(1):I(2),I(3):I(4)),ones(1000,20),BpTau(I(1):I(2),I(3):I(4))];

BW = bwareaopen(BAbeta, 50);

IMG3 = [BW(I(1):I(2),I(3):I(4)),ones(1000,20),BMAP2(I(1):I(2),I(3):I(4)),ones(1000,20),BGFAP(I(1):I(2),I(3):I(4)),ones(1000,20),BLIP(I(1):I(2),I(3):I(4)),ones(1000,20),BDAPI(I(1):I(2),I(3):I(4)),ones(1000,20),BpTau(I(1):I(2),I(3):I(4))];

imwrite([IMG;ones(20,size(IMG,2));IMG1;ones(20,size(IMG,2));IMG2;ones(20,size(IMG,2));IMG3],[filename(1:end-4),'_test.png'])

imwrite(BW,[filename(1:end-4),'_Abeta.png'])
imwrite(BpTau,[filename(1:end-4),'_pTau.png'])

Segmentations.Abeta = BW;
Segmentations.MAP2 = BMAP2;
Segmentations.GFAP = BGFAP;
Segmentations.Lipofuscin = BLIP;
Segmentations.DAPI = BDAPI;
Segmentations.pTau = BpTau;

save([filename(1:end-4),'_segmentation.mat'], '-struct', 'Segmentations')		 

end