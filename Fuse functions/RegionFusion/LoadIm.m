function Im = LoadIm(FileName)
%% Loads image stored in FileName and stores in Im.
%% Converts RGB to grey scale
%% Converts from indexed image to gray scale


[Im,ma] = imread(FileName);

if isind(Im) & ~isempty(ma)  
    Im = double(ind2gray(Im,ma));
else
    if isgray(Im)
        Im = double(Im);
    else
        Im = double(rgb2gray(Im));
    end; 	
end;	