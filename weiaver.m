function weiaver(img1,img2,nr)

im1=double(imread(img1));
im2=double(imread(img2));

file = img1(1:(find(img1 == '.')-1));
ext = img1((find(img1 == '.')+1):end);

weight = (1:nr)/(nr+1)

for i=1:nr
    imw= weight(i).*im1 + (1-weight(i)).*im2;
    imwrite(uint8(imw),[file '_ib_' num2str(i) '.' ext]);
end