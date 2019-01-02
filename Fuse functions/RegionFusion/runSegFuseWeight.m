%%Performing region fusion
clear all;
%Set vars:
Levels = 4;
seg_type = 'joint'; %or 'unimode'
Priority_Method = 1;
Fusion_Method = 1;
Low_Fusion_Method = 0;
Weight = [1, 51, 1.0; 2, 33, -1.5; 1, 49, 1.8];

%Load Ims
Im(:,:,1) = LoadIm('d:\images\2d\mmsegfus\clocks\clockA.tif');
Im(:,:,2) = LoadIm('d:\images\2d\mmsegfus\clocks\clockB.tif');



[overlay,intolay,map,intmap] = Segment_Rob(Im, seg_type, Levels);
imview(overlay(:,:,1),[]);

RegMap = map;   

Fused = Region_Fused_Weight(Im, RegMap, Levels, Priority_Method, Fusion_Method, Low_Fusion_Method, Weight);
imview(Fused,[]);