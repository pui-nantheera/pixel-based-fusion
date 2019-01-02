 function b=imgaus(a,p3,p4)
if ~isa(a, 'double')
    cl = class(a);
    a=im2double(a);
    chc = 1;
else
    chc = 0;
end
b = a + sqrt(p4)*randn(size(a)) + p3;
 