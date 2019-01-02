function [res_str,titlestr,r]=metric(ims,imout,v,param);
%metrics.m, v 0.1 2004/02/17 AL: ,r
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%     function [results,titlestr]=metric(ims,imout,v,params);
%
%  Inputs:      ims - Input Images to be fused
%             imout - Output Images to be fused
%                 v - Metric (see below)
%             param - Cell Array with Input Parameters (see below)
%                 First column contains the names of the parameters
%                 Second column contains the values of the parameters
%
% Outputs:  results - String with results
%          titlestr - String with information about the performed process
%
%===========================================================================
%
%    Metrics:
%
%      v =  1: Mean of each input image (no input arguments)
%
%      v =  2: Petrovic's Metric (no input arguments)  (petmetric.m for details)
%
%      v =  3: Piella's Metric                         (imqind.m for details)
%                Input parameters needed: bs, cedge
%

% read parameters
for i = 1:size(param,1)
  command = param{i,1};
  command = [command '= (param{i,2});'];
  eval(command);
end

% Remove the comment from the next line to gain access to the command
% window. Then, type params if you want to see the parameters you are
% sending to this function

% keyboard

frames=size(ims,4);
res_str{1,1}=[];
for frame=1:frames
  res_str{end,1}=['frame' num2str(frame)];

  imi=ims(:,:,:,frame,:);
  imo=imout(:,:,:,frame,:);
  % AL:
  if size(imi,3)==3 % ((v <= 4) | (v==6)) &      % these metrics only work with gray images
    imi = rgb2grayn(imi);
    imo = rgb2grayn(imo);
  end

  imi = squeeze(imi);
  imo = squeeze(imo);

  switch v
    case 1, % this is Metric 1 it calculates the mean of input images
      titlestr='METR MEA';
      r=mymetric1(imi(:,:,sele));
      for i=1:size(r,3) % convert r into a cell array of strings
        res_str{end+1,1}=[' im' num2str(sele(i)) ' = ' num2str(r(:,:,i))];
      end
    case 2, % Mutual Information
      titlestr='METR MUT';
      r=mif(imi(:,:,sele),imo);
      res_str{end+1,1}=[' M.Inf = ' num2str(r)];
    case 3, % Petrovic and Xydeas Metric
      titlestr='METR PET';
      [q, q_map]=petmetric(imi(:,:,sele), imo);
      % AL:
      r = q;
      res_str{end+1,1}=[' Petmetric' num2str(frame) ' = ' num2str(q)];
    case 4, % Piella's Quality Index
      titlestr=['METR PIE bs= ' num2str(bs) '; cedge=' num2str(cedge) ';'];
      % AL:
      [q, q_map] = imqmet(imi(:,:,sele), imo, bs, cedge);
      % AL:
      r = q;
      res_str{end+1,1}=[' Qindex' num2str(frame) ' = ' num2str(q)];
    case 5, % AL: DWT SNR
      titlestr=['METR DSN lvl= ' num2str(lvl)];% '; wavf= ' wavf '; recon=0;'];
      r = DWTSNR_loop(imi(:,:,sele),imo,lvl);
      r = mean(r);% ??
      res_str{end+1,1}='DWT SNR';
      for i=1:size(r,3) % convert r into a cell array of strings
        res_str{end+1,1}=[' im' num2str(sele(i)) ' = ' num2str(r(:,:,i))];
      end
    case 6, % AL: DWT MI
      titlestr='METR DMI'; % titlestr='METR DSF'; % r = DWTSNR_loop(mean(imi(:,:,sele),3),imo); % r = 1./r;
%      r = SubbandMI_DWT(imi(:,:,sele),imo);
      r = QI_DWT(imi(:,:,sele),imo);
      r = r(:,1:end-1); % not using approximations (negligible for large J)
      r = mean(r(:));
      res_str{end+1,1}=[' M.Inf DWT = ' num2str(r)];
    case 7, % AL: SSIM
      titlestr='SSIM'; 
      r = mean( ...
        [ssim_index(255*imi(:,:,sele(1)), 255*imo) ssim_index(255*imi(:,:,sele(2)), 255*imo)] ...
        );
      res_str{end+1,1}=[' SSIM = ' num2str(r)];
    case 8, % AL: CWSSIM
      titlestr='CWSSIM'; 
      r = mean( ...
        [cwssim_DWT(255*imi(:,:,sele(1)), 255*imo) cwssim_DWT(255*imi(:,:,sele(2)), 255*imo)] ...
        );
      res_str{end+1,1}=[' CWSSIM = ' num2str(r)];
    case 9, % AL: PSNR
      titlestr='PSNR'; 
       r = (psnr(imo,imi(:,:,sele(1))) + psnr(imo,imi(:,:,sele(2))))/2;
       res_str{end+1,1}=[' PSNR = ' num2str(r)];
    otherwise,
  end
  res_str{end+1,1}=[];
end

titlestr = [titlestr ' sele=' num2str(sele) ';'];