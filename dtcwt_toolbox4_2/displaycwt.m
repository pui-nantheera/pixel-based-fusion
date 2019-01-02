function displaycwt(filename) 
 
 
m = gray(256);
[M1,ma] = imread(filename);
if isind(M1) & ~isempty(ma)  
    M1 = 256*double(ind2gray(M1,ma));
else
    if isgray(M1)
       M1 = double(M1);
    else
       M1 = double(rgb2gray(M1));
    end; 	
end;	

zt = 4; 
E = cell(6,zt); 
 
 
[C,S] = dtwavedec2(M1,zt+1,'antonini','qshift_c'); 

for i1 = 1:zt 
  % calculate and store actual image size 
 % Pick out 6 subbands at each level and display them. 
	b1 = cwtband6(C,S,i1);    

% select coefficients and store them 
	for i = 1:6  
   	E(i,i1) = {b1(:,:,i)}; 
	end; 
	%{selc(b1, b2, ap);   
   
end; 
 
l4 = cwtband2(C,S,zt+1,'l','real'); 
prev1 = l4; 
prev2 = l4; 
 
for i1 = zt:-1:1 
    
   prev1 = [E{3,i1} E{1,i1} ; E{5,i1} prev1];  
   prev2 = [E{6,i1} prev2 ; E{4,i1} E{2,i1}];
 
end; 
 
outIm = [prev1; prev2]; 
% image(2* abs(outIm)); 
image(4 * abs(outIm)); 
axis image; 
colormap gray; 
figure; 
image(5 * angle(outIm)); 
axis image; 
colormap gray; 
 
