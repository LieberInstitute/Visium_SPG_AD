
function rnascope_human(filename)

img = load(filename);
O = fieldnames(img);
[Y,~,Z] = size(img.(O{1}));

for i = 1:numel(O)
	  %channe_i = ['im2double(img.',O{i},')']; %for rosehip
	  channe_i = ['rescale(img.',O{i},')'];
       
  if contains(channe_i, 'Lipofuscin')
      channel = eval(channe_i);
  elseif contains(channe_i,'DAPI')
	  channel = medfilt3(eval(channe_i),[7 7 3]);
  else
      channel = imhmin(eval(channe_i),std2(eval(channe_i)));% suppress background noise in RNA scope channels.
  end
  
   thresh = graythresh(channel); %for rosehip
   BWc = imbinarize(channel,thresh);
	  
	 %if thresh<0.04
	 %  BWc = imbinarize(channel,0.04);
	 %end
  
if contains(channe_i,'DAPI')
    BWc = imfill(BWc,'holes');
	bw3=max(BWc,[],3);
    D = -bwdist(~bw3);
    mask = imextendedmin(D,2);
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    bw3(Ld2 == 0) = 0;

    setenv('TZ','America/New_York')
parfor zi =1:Z	
		A = BWc(:,:,zi);
		A(bw3==0)=0;
		BWc(:,:,zi) = A;	 	
end
	
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
for pp = 1:numel(O)
IMG = [IMG,[max(mat2gray(img.(O{pp})),[],3),ones(Y,20); max(Segmentations.(O{pp}),[],3), ones(Y,20)]];
end

imwrite(IMG,[filename(1:end-4),'.png']);
end
