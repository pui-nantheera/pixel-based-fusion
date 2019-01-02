function [op,outs]=interpreter(history,op,outs)
%interpreter.m, v 0.1 2004/01/29
%===========================================================================
%               Eduardo Fernandez Canga - University of Bristol
%
%                        Copyright (c) 2004
%===========================================================================
%  function [op,outs] = interpreter(history,op,outs)
%   
%     op   == operation active at the moment
%                     1: Registration
%                     2: Processing
%                     3: Fusion
%                     4: Metrics
%     outs == # of outputs of the operation
%                     0: only one output
%                     1: an output per input
%     history == line of history read
global ims imout

if size(history)<4
    return
end
order=history(1:4);

switch order
    case 'LOAD'
        M1=LoadIm(history(6:end));
        
        if prod(size(ims)) == 1
            ims=M1;
        elseif (size(M1) ~= size(ims(:,:,:,1))) & any(ims(:))
            error('Input Images must be same size','FUSION ERROR','CREATEMODE')
            return
        else
            ims(:,:,:,:,size(ims,5)+1)=M1;
        end     
        
    case 'LOAV'
        y=avi2vol(history(15:end));
        
        if prod(size(ims)) == 1
            ims=y;
        elseif any(size(y) ~= size(ims(:,:,:,:,1))) & any(ims(:))
            error('Input Videos must be same size','FUSION ERROR','CREATEMODE')
            return
        else
            ims(:,:,:,:,size(ims,5)+1)=y;
        end     
        
        %    case 'LOAS' %still to implement
        
    case 'DELE'
        v=str2num(history(12:end));
        c=size(ims,5);
        if c<2 | v==0
            ims=0;
            imout=0;
        else
            r=1:c;
            r=r(r~=v);
            ims=ims(:,:,:,:,r);
        end
        
    case 'UPDA'
        ims=imout;
        imout=0;
        
    case 'ADDB'
        p=size(ims,5)+(1:size(imout,5));
        ims(:,:,:,:,p)=imout;
        
    case 'HIST' %do nothing 
        
    case 'DIST'
        op   = 2;
        outs = 1;
        
        meth=history(6:9);
        v=get_dist_meth(meth);
        param = text2cell(history(10:end));
        
        imout = proc(ims,v,param);
        
    case 'FUSE'
        op   = 3;
        outs = 0;
      
        meth=history(6:8);
        v=get_fuse_meth(meth);
        param = text2cell(history(9:end));

        [imout,a]=fuse(ims,v,param);        
                
    case 'METR'
        op = 4;

        meth=history(6:8);
        v=get_metr_meth(meth);
        param = text2cell(history(9:end));
        
        [res_str]=metric(ims,imout,v,param);
        disp(res_str)
        % How to pass this output to the interface?? New window??
    case 'COLO'
        op = 2;
        outs =1;
        param = text2cell(history(5:end));
        
        imout = proc(ims,31,param);
        
        
        
       
    otherwise
        txt=sprintf(['\nCan''t read ' order ' command:\n\n%s'], history);
        error(txt)
        return
end


function val=get_fuse_meth(v)
switch v
case 'AVR', val = 1;
case 'MAX', val = 3;
case 'MIN', val = 4;
case 'CEM', val = 5;
case 'LAP', val = 6;
case 'FSD', val = 7;
case 'RAT', val = 8;
case 'CON', val = 9;
case 'GRA', val = 10;    
case 'DWB', val = 11;    
case 'SIH', val = 12;    
case 'MOD', val = 13;    
case 'CWT', val = 14;        
case 'LIS', val = 15;    
case 'MFO', val = 17;    
case 'REG', val = 18;    
end

function val=get_metr_meth(v)
switch v
case 'MEA', val = 1;
case 'MUT', val = 2;
case 'PET', val = 3;
case 'PIE', val = 4;
case 'DSN', val = 5;
end

function val=get_dist_meth(v)
switch v
case 'GAUS', val = 11;
case 'SALT', val = 12;
case 'SPEC', val = 13;
case 'POIS', val = 14;
end


function param = text2cell(strs)
% the first character passed must be a blank character
p = find(strs=='=');
q = [0 find(strs==';')];

for i=1:size(p,2)
    param{i,1}=strs((q(i)+2):(p(i)-1));
    param{i,2}=[str2num(strs((p(i)+1):(q(i+1)-1)))];
    if isempty(param{i,2})
        param{i,2}=strs((p(i)+1):(q(i+1)-1));
    end
end
