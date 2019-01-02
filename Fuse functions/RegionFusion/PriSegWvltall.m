function Priority = PriSegWvltall(HighCoef, LabDecomp, Levels, ImSize, Priority_Method)
%% This function calcultes the priority over the image depending on various
%% factors such as entropy. Priority is the same for a region over all 
%% subbands and levels.  
%% HighCoef should be a cell array containing Levels number of cells.  Each
%% cell should contain the high pass coefficients for the 6 complex wavelet
%% orientations.  
%% LabDecomp contains a labeled map of the segmentation of the image.
%% Priority_method is the type of priority calculated where:   
%%                  1 - Entropy based priority
%%                  2 - Variance based
%%                  3 - Activity based


lab = 1;
Priority = cell(Levels, 1);
levsize = ceil(ImSize./2);
for a = 1:Levels
    Priority(a, :) = { zeros(levsize(1), levsize(2), 6) };    
    levsize = ceil(levsize./2);
end;

while (any(any(LabDecomp{1}==lab)))
    levsize = ceil(ImSize./2);
    prient = 0; 
    privar = 0;
    priact = 0;
    for a = 1:Levels
        seg = find(LabDecomp{a} == lab);
        if (length(seg) > 0) 
            for b = 1:6
                switch Priority_Method
                    case 1
                        %% Calculate the entropy
                        temp = HighCoef{a}(:,:, b);
                        temp1 = abs(temp(seg)).^2; 
                        temp2 = temp1.*log(temp1);
                        temp2(isnan(temp2)) = 0;
                        prient = prient + ((sum(sum(temp2)))/size(seg,1));
                                
                        clear temp1;
                        clear temp2;
                    case 2              
                        %% Calculate variance based priority
                        temp = HighCoef{a}(:,:, b);
                        privar = privar + var(temp(seg), 1);
                    case 3
                        %% Calculate activity based priority
                        temp = HighCoef{a}(:,:, b);
                        priact = priact + sum(sum(abs(temp(seg))))/size(seg, 1);
                        
                end;
            end;
            clear temp;
            clear seg;
        end;
        levsize = ceil(levsize./2);        
    end;
    
    %assign values
    for a = 1:Levels
        seg = find(LabDecomp{a} == lab);
        if (length(seg) > 0) 
            for b = 1:6
                switch Priority_Method
                    case 1                
                        temp3 = Priority{a}(:,:, b);
                        temp3(seg) = prient/(Levels*6);
                        Priority{a}(:,:, b) = temp3;
                    case 2
                        temp3 = Priority{a}(:,:, b);
                        temp3(seg) = privar/(Levels*6);
                        Priority{a}(:,:,b) = temp3;
                    case 3
                        temp3 = Priority{a}(:,:, b);
                        temp3(seg) = priact/(Levels*6);
                        Priority{a}(:,:,b) = temp3;
                end;
            end;
            clear temp3;
        end;
        levsize = ceil(levsize./2);        
    end;
    
    lab = lab + 1;
end;
