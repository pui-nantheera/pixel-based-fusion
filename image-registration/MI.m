
function h=MI(image_1,image_2)

a=joint_h(image_1,image_2);
[r,c] = size(a);
b= a./(r*c);
y_marg=sum(b);
x_marg=sum(b');

Hy=0;
for i=1:c;    %  col
      if( y_marg(i)==0 )
         %do nothing
      else
         Hy = Hy + -(y_marg(i)*(log2(y_marg(i))));
      end
   end
   
Hx=0;
for i=1:r;    %rows
   if( x_marg(i)==0 )
         %do nothing
      else
         Hx = Hx + -(x_marg(i)*(log2(x_marg(i))));
      end   
   end
h_xy = -sum(sum(b.*(log2(b+(b==0)))));
h = Hx + Hy - h_xy;




