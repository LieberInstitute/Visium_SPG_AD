toolbox='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/code/VistoSeg/';
addpath(genpath(toolbox))

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD1_V10A27-004_A1.mat';
%segmentIF(filename,[8301,9300,9501,10500])
refineIF(filename,[8301,9300,9501,10500],[],[],0.2,[])

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD1_V10A27-004_D1.mat';
%segmentIF(filename,[4001,5000,11501,12500])
refineIF(filename,[4001,5000,11501,12500],0.04,0.2,0.2,0.2)
	
filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_A1.mat';
%segmentIF(filename,[12501,13500,5801,6800])
refineIF(filename,[12501,13500,5801,6800],[],[],0.2,[])

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_B1.mat';
%segmentIF(filename,[8701,9700,6801,7800])
refineIF(filename,[8701,9700,6801,7800],0.2,0.2,0.25,0.2)

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_C1.mat';
%segmentIF(filename,[15001,16000,7601,8600])
refineIF(filename,[15001,16000,7601,8600],0.06,0.2,0.25,[])

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD2_V10A27-106_D1.mat';
%segmentIF(filename,[4801,5800,8001,9000])
refineIF(filename,[4801,5800,8001,9000],0.15,0.3,0.3,0.2)

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_A1.mat';
%segmentIF(filename,[13501,14500,8701,9700])
refineIF(filename,[13501,14500,8701,9700],[],[],0.15,0.1)

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_B1.mat';
%segmentIF(filename,[7001,8000,9501,10500])
refineIF(filename,[7001,8000,9501,10500],0.1,[],0.15,[])

filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_C1.mat';
%segmentIF(filename,[13001,14000,10201,11200])
refineIF(filename,[13001,14000,10201,11200],0.1,0.25,0.3,0.2)


filename='/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/VIFAD3_V10T31-036_D1.mat';
%segmentIF(filename,[4251,5250,7501,8500])
refineIF(filename,[4251,5250,7501,8500],0.2,0.3,0.4,0.25)

movefile('/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/*.png','/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Segmentations/')
movefile('/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Capture_Areas/*segmentation.mat','/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/Images/VistoSeg/Segmentations/')

