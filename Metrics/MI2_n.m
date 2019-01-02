function h=MI2(image_1,image_2)
% function h=MI2(image_1,image_2)
%
% Takes a pair of images and returns the mutual information Ixy using joint entropy function JOINT_H.m
% 
% written by http://www.flash.net/~strider2/matlab.htm

a=joint_h(image_1,image_2); % calculating joint histogram for two images
[r,c] = size(a);
b= a./(2*r*c); % normalized joint histogram
y_marg=sum(b); %sum of the rows of normalized joint histogram
x_marg=sum(b');%sum of columns of normalized joint histogran

% x_marg=x_marg(ones(1,256),:);
% y_marg=y_marg(ones(1,256),:);

% h=0;
% for i=1:256
%     for j=1:256
%         if (x_marg(j)*y_marg(i)) ~= 0 & (b(i,j) ~=0)
%             h=h+(b(i,j)*log(b(i,j)/(x_marg(j)*y_marg(i))));
%         end
%     end
% end


Hy=0;
for i=1:c;    %  col
     if( y_marg(i)==0 )
        %do nothing
     else
        Hy = Hy + -(y_marg(i)*(log(y_marg(i)))); %marginal entropy for image 1
     end
   end
   
Hx=0;
for i=1:r;    %rows
   if( x_marg(i)==0 )
         %do nothing
      else
         Hx = Hx + -(x_marg(i)*(log(x_marg(i)))); %marginal entropy for image 2
      end   
   end
% Hx = -sum(sum(b.*(log(x_marg+(x_marg==0)))));
% Hy = -sum(sum(b.*(log(y_marg+(y_marg==0)))));
h_xy = -sum(sum(b.*(log(b+(b==0))))); % joint entropy
h = Hx + Hy - h_xy;% Mutual information