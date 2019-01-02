function met=metr_compr_im(image_file,image_file2,compress)

block_size=8;
j=1;
for i=compress
    image_file1=[image_file(1:find(image_file=='.')-1) '_' sprintf('%.3d',i) '.jpg'];
    im=double(imread(image_file));
    im(:,:,2)=double(imread(image_file2));
    
    imf=double(imread(['f' image_file1]));
    
    [met(j),map]=imqind(im, imf, block_size);
    met(j)
    j=j+1;
end


