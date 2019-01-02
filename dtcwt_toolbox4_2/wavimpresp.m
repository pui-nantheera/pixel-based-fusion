function outvol = wavimpresp(odd_filt,q_filt,level);

outvol =cell(level,1);

load(odd_filt);
load(q_filt);

H0t = [h0o;0] + j*[0;h0o];
H1t = [h1o;0] + j*[0;h1o];

%pad everything to greater length
Hlength = max(length(H0t),length(H1t));
H0 = zeros(Hlength,1);
H1 = H0;

%choose insertion points to centre
insert1 = fix((Hlength-length(H1t))./2) +1;
insert0 = fix((Hlength-length(H0t))./2) +1;
%insert into padding
H1(insert1:insert1+length(H1t)-1) = H1t;
H0(insert0:insert0+length(H0t)-1) = H0t;
%           column      row
%               |         |
%               v         v
LL(:,:,1) = real(H0)*real(H0).';
LL(:,:,2) = real(H0)*imag(H0).';
LL(:,:,3) = imag(H0)*real(H0).';
LL(:,:,4) = imag(H0)*imag(H0).';

HL(:,:,1) = real(H1)*real(H0).';
HL(:,:,2) = real(H1)*imag(H0).';
HL(:,:,3) = imag(H1)*real(H0).';
HL(:,:,4) = imag(H1)*imag(H0).';

LH(:,:,1) = real(H0)*real(H1).';
LH(:,:,2) = real(H0)*imag(H1).';
LH(:,:,3) = imag(H0)*real(H1).';
LH(:,:,4) = imag(H0)*imag(H1).';

HH(:,:,1) = real(H1)*real(H1).';
HH(:,:,2) = real(H1)*imag(H1).';
HH(:,:,3) = imag(H1)*real(H1).';
HH(:,:,4) = imag(H1)*imag(H1).';

%gen first level highpass outputs
j2 = sqrt([0.5 -0.5]);
outvol{1}(:,:,1) = (HL(:,:,1) - HL(:,:,4)).*j2(1) + (HL(:,:,2) + HL(:,:,3)).*j2(2);
outvol{1}(:,:,2) = (HL(:,:,1) + HL(:,:,4)).*j2(1) + (HL(:,:,2) - HL(:,:,3)).*j2(2);
outvol{1}(:,:,3) = (LH(:,:,1) - LH(:,:,4)).*j2(1) + (LH(:,:,2) + LH(:,:,3)).*j2(2);
outvol{1}(:,:,4) = (LH(:,:,1) + LH(:,:,4)).*j2(1) + (LH(:,:,2) - LH(:,:,3)).*j2(2);
outvol{1}(:,:,5) = (HH(:,:,1) - HH(:,:,4)).*j2(1) + (HH(:,:,2) + HH(:,:,3)).*j2(2);
outvol{1}(:,:,6) = (HH(:,:,1) + HH(:,:,4)).*j2(1) + (HH(:,:,2) - HH(:,:,3)).*j2(2);

%gen second level filters, remembering to zero-insert
H00a = dyadup(h0a,0);
H01a = dyadup(h1a,0);
H00b = dyadup(h0b,0);
H01b = dyadup(h1b,0);

%do further level filtering
for n=2:level;
    %LL - col then row
    clear newLL;
    temp = conv2(LL(:,:,1),H00a,'full');
    newLL(:,:,1) = conv2(temp,H00a.','full');
    temp = conv2(LL(:,:,2),H00a,'full');
    newLL(:,:,2) = conv2(temp,H00b.','full');
    temp = conv2(LL(:,:,3),H00b,'full');
    newLL(:,:,3) = conv2(temp,H00a.','full');
    temp = conv2(LL(:,:,4),H00b,'full');
    newLL(:,:,4) = conv2(temp,H00b.','full');
    
    clear LH HL HH;
    
    %HL - col then row
    temp = conv2(LL(:,:,1),H01a,'full');
    HL(:,:,3) = conv2(temp,H00a.','full');
    temp = conv2(LL(:,:,2),H01a,'full');
    HL(:,:,4) = conv2(temp,H00b.','full');
    temp = conv2(LL(:,:,3),H01b,'full');
    HL(:,:,1) = conv2(temp,H00a.','full');
    temp = conv2(LL(:,:,4),H01b,'full');
    HL(:,:,2) = conv2(temp,H00b.','full');
    
    %LH - col then row
    temp = conv2(LL(:,:,1),H00a,'full');
    LH(:,:,2) = conv2(temp,H01a.','full');
    temp = conv2(LL(:,:,2),H00a,'full');
    LH(:,:,1) = conv2(temp,H01b.','full');
    temp = conv2(LL(:,:,3),H00b,'full');
    LH(:,:,4) = conv2(temp,H01a.','full');
    temp = conv2(LL(:,:,4),H00b,'full');
    LH(:,:,3) = conv2(temp,H01b.','full');
    
    %HH - col then row
    temp = conv2(LL(:,:,1),H01a,'full');
    HH(:,:,4) = conv2(temp,H01a.','full');
    temp = conv2(LL(:,:,2),H01a,'full');
    HH(:,:,3) = conv2(temp,H01b.','full');
    temp = conv2(LL(:,:,3),H01b,'full');
    HH(:,:,2) = conv2(temp,H01a.','full');
    temp = conv2(LL(:,:,4),H01b,'full');
    HH(:,:,1) = conv2(temp,H01b.','full');
    
    LL = newLL;
    
    outvol{n}(:,:,1) = (HL(:,:,1) - HL(:,:,4)).*j2(1) + (HL(:,:,2) + HL(:,:,3)).*j2(2);
    outvol{n}(:,:,2) = (HL(:,:,1) + HL(:,:,4)).*j2(1) + (HL(:,:,2) - HL(:,:,3)).*j2(2);
    outvol{n}(:,:,3) = (LH(:,:,1) - LH(:,:,4)).*j2(1) + (LH(:,:,2) + LH(:,:,3)).*j2(2);
    outvol{n}(:,:,4) = (LH(:,:,1) + LH(:,:,4)).*j2(1) + (LH(:,:,2) - LH(:,:,3)).*j2(2);
    outvol{n}(:,:,5) = (HH(:,:,1) - HH(:,:,4)).*j2(1) + (HH(:,:,2) + HH(:,:,3)).*j2(2);
    outvol{n}(:,:,6) = (HH(:,:,1) + HH(:,:,4)).*j2(1) + (HH(:,:,2) - HH(:,:,3)).*j2(2);
    
    H00a = dyadup(H00a,0);
    H01a = dyadup(H01a,0);
    H00b = dyadup(H00b,0);
    H01b = dyadup(H01b,0);
end;