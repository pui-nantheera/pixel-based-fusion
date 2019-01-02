function Fused = fuse_reg(Im, RegMap, Levels, Priority_Method, Fusion_Method, Low_Fusion_Method, Weight)

%% This function takes N registered greyscale images, Im, with either a single
%% joint segmentation map or seperate segmentation maps, RegMaps, and fuses them 
%% depending on the remaining parameters.
%%
%% Priority Method is the choice of method of generating priority
%%          1 - Entropy
%%          2 - Variance
%%          3 - Average Activity
%%
%% Fusion_Method defines how the CWT coef. are combined
%%          1 - Take all coefficients from region with highest priority
%%          2 - Take a average of coefficients weighted according to
%%              priority
%%
%% Low_Fusion_Method defines how low pass coefficients are combined
%%          0 - Take average
%%          1 - Take all values from Image 1.
%%          2 - Take all values from Image 2. etc.
%%
%% Weight is a matrix containing [Im1, reg1, weight1; Im2,reg2, weight2;.....]
%%          Where Im_ is the image number (as in the Im matrix) conatining
%%          the region, reg_ to be weighted with weighting weight_;
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                    %%  
%%          by John Lewis                                                             %%
%%          University Of Bristol, UK                                                 %%
%%          Copyright 2003                                                            %%
%%                                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Weight(1,1:3) == [-1,-1, -1]
    Weight = [1 1 1];
end

ImSize = size(Im);

%% Check all input parameters valid
if ((size(RegMap,1) ~= ImSize(1)) | (size(RegMap,2) ~= ImSize(2)))
    error('Region map must be the same size as the image');
end;

%% expany joint segmentation map so that it is copied to make a volume the
%% same size as the image
if  (size(RegMap,3) ~= ImSize(3))
    if (size(RegMap,3) == 1)
        for i = 2:ImSize(3)
            RegMap(:,:,i) =  RegMap(:,:,1);
            warning('Joint segmentation assumed and used for each input image'); 
        end;
    else
        error('Region map must be a volume the same size as the image volume or a single (joint) segmentation');        
    end;
end;

if ((Priority_Method ~= 1) & (Priority_Method ~= 2) & (Priority_Method ~= 3))
    error('Invalid Priority_Method');
end;

if ((Fusion_Method ~= 1) & (Fusion_Method ~= 2))
    error('Invalid Fusion Method');
end;

if ((Low_Fusion_Method < 0) | (Low_Fusion_Method > ImSize(3)))
    error('invalid Low_Fusion_Method')
end;

if ( size(Weight, 2) ~= 3 )
    error('weight must have 3 columns');    
end;

for i = 1:size(Weight,1)
    if (Weight(i,1) < 1) | (Weight(i,1) > ImSize(3)) 
        error('invalid image numbers in Weight. Should have the form [Im1, reg1, Weight1; Im2,reg2, Weight2;.....]');
    end;
    if (Weight(i,2) < 1) | (Weight(i,2) > max(max(RegMap(:,:,Weight(i,1)))) ) 
        error('invalid region number in weight. Should have the form [Im1, reg1, weight1; Im2,reg2, weight2;.....]');
    end;
end;



Priority = cell(Levels, ImSize(3));
highcoef = cell(Levels, ImSize(3));


%Cycle through all images calculating data
for a = 1:ImSize(3)
    %% Calculate the Dual-Tree Complex Wavelet transform coefficients
    [lowcoef(:,:,a),highcoef(:,a)] = dtwavexfm2(Im(:,:,a),Levels,'antonini','qshift_06');

    %find the corresponding label map at each level 
    RegDecomp(:,a) = LabDecompProt(RegMap(:,:,a), Levels, ImSize);
    %% calculate priority maps for each subband at each level
    %% Priority returns 3 different priority maps, Priority{level, type}
    %% where type is:   1 - Entropy based priority
    %%                  2 - Variance based
    %%                  3 - Activity based
    Priority(:,a) = PriSegWvltall(highcoef(:,a), RegDecomp(:,a), Levels, ImSize, Priority_Method);

    %% If average low pass methods
    if (a==1)
        Fused_LowCoef = zeros(size(lowcoef(:,:,1)));
    end;
    if (Low_Fusion_Method == 0)
        Fused_LowCoef = Fused_LowCoef + lowcoef(:,:,a);
    end;
end;

%% Fuse Low pass Data
if (Low_Fusion_Method == 0)
    %% just average the low pass coefficients together for low pass
    %% coefficeints    
    Fused_LowCoef = Fused_LowCoef./ImSize(3);
else
    %% Otherwise takes vales from appropriate image
    Fused_LowCoef = lowcoef(:,:,Low_Fusion_Method);
end;

%add a weighting eg. to the brighter parts of IR image
WeightMask = cell(Levels, ImSize(3));
levsize = [ceil(ImSize(1)./2) ceil(ImSize(2)./2)]; 

for a = 1:Levels
    WeightMask(a, :) = { zeros(levsize(1), levsize(2)) };    
    levsize = ceil(levsize./2);
end;

Weight(:,3) = Weight(:,3) - 1;
for i = 1:size(Weight,1)
    for j = 1:Levels
        WeightMask{j,Weight(i,1)} = WeightMask{j,Weight(i,1)} + ( RegDecomp{j,Weight(i,1)} == Weight(i,2) ).*Weight(i, 3);
    end;
end;



%% choose the best match of priority at each level at each subband
Mask = cell(Levels,1);
Fused_HighCoef = cell(Levels, 1);

for a = 1:Levels
    Fused_HighCoef{a} = zeros(size(highcoef{a}));

    %add the weighting to all the coefficients
    for i = 1:ImSize(3)
        WeightMask{a,i} = WeightMask{a,i} + 1;
        for j = 1:6
            highcoef{a,i}(:,:,j) = highcoef{a,i}(:,:,j).*WeightMask{a,i};
        end;
    end;

    %Convert from cell to array for max argument
    for b = 1:ImSize(3)
        Prioritytemp1(:,:,:,b) = Priority{a, b};
        HighCoeftemp(:,:,:,b) = highcoef{a, b};
    end;

   
    switch Fusion_Method
        case 1  %% Choose which image regions should come from based on highest priority
            [C Mask{a}] = max(Prioritytemp1,[],4);
            %% Need to add code for when have the same priorities
            for b = 1:ImSize(3)
                Fused_HighCoef{a} = Fused_HighCoef{a} + double(highcoef{a,b}.*(Mask{a}==b));%+highcoef{a,2}.*(Mask{a, 1}==2)%+(max(highcoef{a,1}, highcoef{a,2})).*(Mask{a, 1}==3);        
            end;
        case 2 %% Use a weighted average of the coefficients in each region, depending on priority of that region
            maxprio = zeros(size(Priority{a,1}));
            for b = 1:ImSize(3);
                maxprio = maxprio + Priority{a,b};
            end;
            
            for b = 1:ImSize(3);
                Mask{a} = Priority{a,b}./maxprio;
                Fused_HighCoef{a} = Fused_HighCoef{a} + Mask{a}.*highcoef{a,b};
                
            end;
            clear maxprio;
            
    end;
    clear Prioritytemp1; clear C; clear HighCoeftemp;
end;

Fused = dtwaveifm2(Fused_LowCoef, Fused_HighCoef,'antonini','qshift_06',ones(6,length(Fused_HighCoef)));


