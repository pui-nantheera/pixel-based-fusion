%X = imread('ll4_SDa.png');
X1 = imread('im1_reg_ll3_SD.png');
X2 = imread('im2_ll3_IR.png');

X1 = im2double(imread('C:\Fusion\codes\IFT update\data\clocks\clockA.jpg'));
X2 = im2double(imread('C:\Fusion\codes\IFT update\data\clocks\clockB.jpg'));


fmeth = 'SAS';
J = 0;
enh_meth = 'VBE';

%enh_meth = 'HEQ'
X1 = rgb2gray(X1);
X1 = double(X1);

X2 = rgb2gray(X2);
X2 = double(X2);
%%[y Nsigfused Nsig] = contr_enh(X,fmeth,J,enh_meth);
[y, Nsigfused, Nsig] = FUSEcontr_enh(X1,X2,fmeth,3, enh_meth);
    