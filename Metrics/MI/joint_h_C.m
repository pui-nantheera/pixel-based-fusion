
function h=joint_h_C(image_1,image_2)

image_1= image_1(:);
image_2= image_2(:);
N=256;

% rows=size(image_1);
% h=zeros(N,N);
% for i=1:rows;
%     h(image_1(i)+1,image_2(i)+1)= h(image_1(i)+1,image_2(i)+1)+1;
% end
% figure, tfrplot(h/max(h(:)))

% AL:
h = hist3([image_1 image_2], [N N]);
% figure, tfrplot(h/max(h(:)))
