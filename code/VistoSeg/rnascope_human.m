
function rnascope_human(filename,I)

img = load(filename);
O = fieldnames(img);

for i = 1:numel(O)
  channel = rescale(img.(O{i}));
       
  if contains((O{i}),'DAPI')
	 channel = medfilt3(channel,[7 7 3]);
     thresh = graythresh(channel);
  elseif contains((O{i}),{'Lipofuscin','Abeta'})
     thresh = graythresh(channel);
  else
     thresh = 0.7;
  end
  
  BWc = imbinarize(channel,thresh);
  
  if contains((O{i}),'DAPI')
    BWc = imfill(BWc,'holes');
	bw3=max(BWc,[],3);
    D = -bwdist(~bw3);
    mask = imextendedmin(D,2);
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw3(Ld2 == 0) = 0;
    BWc(bw3==0)=0;	
  else	 	 
	 x = imcomplement(channel);
	 x = imhmin(x,2*std(channel(:)));
	 L = watershed(x);
	 BWc(L==0) = 0;  
  end	 

  [~,no_of_dots] = bwlabeln(BWc);
  disp(['segmented ',O{i}, ': ',num2str(no_of_dots)])
  v = ['Segmentations.',O{i}];
  eval([v '= BWc;']);

end

save([filename(1:end-4),'_segmentation.mat'], '-struct', 'Segmentations')		 

O = fieldnames(img);
IMG = [];
%I = [4001,5000,7101,8100];
for pp = 1:numel(O)
    img1 = rescale(img.(O{pp}));
    IMG = [IMG,[mat2gray(img1(I(1):I(2),I(3):I(4))),ones(1000,20); Segmentations.(O{pp})(I(1):I(2),I(3):I(4)), ones(1000,20)]];
end

imwrite(IMG,[filename(1:end-4),'.png']);
end
