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

for nx=1:nbinsx
    for ny=1:nbinsy
        if h(nx,ny)~=0
            logf = h(nx,ny) * log2( h(nx,ny) / (hx(nx) * hy(ny )));
        else
            logf = 0;
        end;
        I=I+logf;
    end;
end;

% H = h .* log2( h + h==0);
% Hx= hx .* log2( hx + hx==0);
% Hx= hy .* log2( hy + hy==0);




count=sum(h(:));

I=I/count;





