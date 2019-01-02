function I=mi(x,y,descriptor)

if nargin==2
 [h,descriptor]=histogram2(x,y);
elseif nargin==3
 nbins=descriptor(:,1);
 lower=descriptor(:,2);
upper=descriptor(:,3);
[h,descriptor]=histogram2(x,y,nbins,lower,upper);
else
 error('Not correct number of input parameters');
end

h=h/sum(h(:));
hy=sum(h);
hx=sum(h');
nbinsx=descriptor(1,3);
nbinsy=descriptor(2,3);

I=0;
count=0;
warn=0;
for nx=1:nbinsx
  for ny=1:nbinsy
    if h(nx,ny)~=0 
      logf=log2(h(nx,ny)/(hx(nx)*hy(ny)));
if logf<0
 warn=1;
end
    else
      logf=0;
    end;
    count=count+h(nx,ny);
    I=I+h(nx,ny)*logf;
  end;
end;

I=I/count;





