%%% DEMO %%%
clear all
path ( genpath(pwd) , path )

% -------------------------------------------------------------------------
% read images
% -------------------------------------------------------------------------
% example: read all images in specific folder
folder = 'C:\Fusion\codes\ift update\data\UN camp\';
files  = dir([folder, '*.bmp']);
totalImg = length(files);
% read image info;
origImgs  = cell(1,totalImg);
for k = 1:totalImg
    origImgs{k} = im2double(imread([folder, files(k).name]));
end

% output directory
outfolder = 'C:\Fusion\codes\results\UN camp\';
mkdir(outfolder);
% -------------------------------------------------------------------------
% image fusion process
% -------------------------------------------------------------------------
methodList = {'Average','PCA',...
                'MAXIMUM','MINIMUM','SpatFreq',...
                'LAPLACIAN','FSD','RATIO','CONTRAST','GRADIENT','Morphological','DWT','LISQ',...
                'SIDWT',...
                'CWT'};
% for measuring performance
ctime            = zeros(1,length(methodList));
quality_mi       = zeros(1,length(methodList));
quality_minorm   = zeros(1,length(methodList));
quality_petrovic = zeros(1,length(methodList));
quality_piella   = zeros(1,length(methodList));
quality_cvejic   = zeros(1,length(methodList));
% parameter
low = 0;        % average low-pass
high = [1 3];   % max magnitude high-pass
blockSize = 8;  % blocksize for Spatial Frequency
threshold = 1;  % threshold for Spatial Frequency
% running
for method = 1:length(methodList)
    
    [fusedImg,~,ctime(method)] = fuse(origImgs, methodList{method},'high',high,'low',low,...
        'blockSize',blockSize,'threshold',threshold);
%     imwrite(fusedImg, [outfolder, 'fusedImg_',methodList{method},'.tif'], 'tif');
%     figure; imshow(fusedImg); title(methodList{method});
    
    % quality assessment
    % ---------------------------------------------------------------------
    % mutual information
    [quality_mi(method),quality_minorm(method)] = mif(origImgs,fusedImg);
    % Petrovic and Xydeas Metric
    quality_petrovic(method) = petmetric(origImgs, fusedImg);
    % Piella's Quality Index
    quality_piella(method) = imqmet(origImgs, fusedImg);
    % Cvejic's Quality Index
    quality_cvejic(method) = Cvejic_metric(origImgs,fusedImg);
end

% -------------------------------------------------------------------------
% model-based image fusion process
% -------------------------------------------------------------------------
modelList = {'WeigtedAvg','FLOM','FLOM-Cauchy','Meridian','Cauchy'};
totalModelMethod = 5;
ctimeM            = zeros(1,totalModelMethod);
qualityM_mi       = zeros(1,totalModelMethod);
qualityM_minorm   = zeros(1,totalModelMethod);
qualityM_petrovic = zeros(1,totalModelMethod);
qualityM_piella   = zeros(1,totalModelMethod);
qualityM_cvejic   = zeros(1,totalModelMethod);
for method = 1:5
    [fusedImg,~,ctimeM(method)] = fuse(origImgs, 'Model','low',low, 'method',method);
%     figure; imshow(fusedImg); title(['Model-based method : ',num2str(method)]);
%     imwrite(fusedImg, [outfolder, 'fusedImg_model',modelList{method},'.tif'], 'tif');
    
    % quality assessment
    % ---------------------------------------------------------------------
    % mutual information
    [qualityM_mi(method),qualityM_minorm(method)] = mif(origImgs,fusedImg);
    % Petrovic and Xydeas Metric
    qualityM_petrovic(method) = petmetric(origImgs, fusedImg);
    % Piella's Quality Index
    qualityM_piella(method) = imqmet(origImgs, fusedImg);
    % Cvejic's Quality Index
    qualityM_cvejic(method) = Cvejic_metric(origImgs,fusedImg);
end
statsList = {'Laplacian','Cauchy','GGaussian','Alpha-Stable'};
for method = 1:4
    [fusedImg,~,ctimeM(totalModelMethod+method)] = fuse(origImgs, 'Stats','low',low, 'method',method);
%     figure; imshow(fusedImg); title(['Model-based method : ',num2str(method)]);
%     imwrite(fusedImg, [outfolder, 'fusedImg_model',modelList{method},'.tif'], 'tif');
    
    % quality assessment
    % ---------------------------------------------------------------------
    % mutual information
    [qualityM_mi(totalModelMethod+method),qualityM_minorm(totalModelMethod+method)] = mif(origImgs,fusedImg);
    % Petrovic and Xydeas Metric
    qualityM_petrovic(totalModelMethod+method) = petmetric(origImgs, fusedImg);
    % Piella's Quality Index
    qualityM_piella(totalModelMethod+method) = imqmet(origImgs, fusedImg);
    % Cvejic's Quality Index
    qualityM_cvejic(totalModelMethod+method) = Cvejic_metric(origImgs,fusedImg);
end

% -------------------------------------------------------------------------
% write objective results
% -------------------------------------------------------------------------
QM = [quality_mi' quality_minorm' quality_petrovic' quality_piella' quality_cvejic'];
QM = [QM; qualityM_mi' qualityM_minorm' qualityM_petrovic' qualityM_piella' qualityM_cvejic'];
CT = [ctime'; ctimeM'];
resultsMat = [(1:24)' QM sum(QM(:,3:5),2) CT];
methodNames = [methodList'; modelList'; statsList'];

fileID = fopen([outfolder,'qualtymetric.txt'],'w');
fprintf(fileID,'%15s %4s %8s %8s %8s %8s %8s %8s %8s\r\n','Method','ID','MI','MInorm', 'Petrovic', 'Pialla', 'Cvejic', 'sumLast3', 'ctime');
for k = 1:length(methodNames)
    fprintf(fileID,'%15s %4u %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f\r\n', char(methodNames{k}), resultsMat(k,:));
end
fclose(fileID);