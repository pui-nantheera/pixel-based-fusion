function b=imsaltpep(a,p3)
if ~isa(a, 'double')
    cl = class(a);
    a=im2double(a);
    chc = 1;
else
    chc = 0;
end
b = a;
x = rand(size(a));
d = find(x < p3/2);
b(d) = 0; % Minimum value
d = find(x >= p3/2 & x < p3);
b(d) = 1;

