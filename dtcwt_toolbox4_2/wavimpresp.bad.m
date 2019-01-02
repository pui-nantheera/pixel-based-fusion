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

LL(:,:,1) = H0*H0.';
LL(:,:,2) = H0*H0'; %row conjugate

if level>=1;
    outvol{1}(:,:,1) = H0*H1.';
    outvol{1}(:,:,2) = H0*H1';
    outvol{1}(:,:,3) = H1*H0.';
    outvol{1}(:,:,4) = H1*H0';
    outvol{1}(:,:,5) = H1*H1.';
    outvol{1}(:,:,6) = H1*H1';
end;

H00 = [h0a] + j*[h0b];
H01 = [h1a] + j*[h1b];

    H00 = dyadup(H00);
    H01 = dyadup(H01);

for n=2:level;
    H00 = dyadup(H00);
    H01 = dyadup(H01);
    collow1 = conv2(LL(:,:,1),H00,'full');
    collow2 = conv2(LL(:,:,2),H00,'full');
    lowlow1 = conv2(collow1,H00','full'); %conj version
    lowlow2 = conv2(collow2,H00.','full');
    colhigh1 = conv2(LL(:,:,1),H01,'full');
    colhigh2 = conv2(LL(:,:,2),H01,'full');
    outvol{n}(:,:,1) = conv2(collow1,H01.','full');
    outvol{n}(:,:,2) = conv2(collow2,H01','full');
    outvol{n}(:,:,3) = conv2(colhigh1,H00.','full');
    outvol{n}(:,:,4) = conv2(colhigh2,H00','full');
    outvol{n}(:,:,5) = conv2(colhigh1,H01.','full');
    outvol{n}(:,:,6) = conv2(colhigh2,H01','full');
    LL = cat(3,lowlow1,lowlow1);
end;

