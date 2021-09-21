function masking(filename)
seg = load(filename);
O = fieldnames(seg);
[Y,~,Z] = size(seg.(O{1}));

LIP = find(contains(O,'Lipofuscin'));  



