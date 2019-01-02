function fftim = wavfreq(invol,padsz);
fftim = zeros(padsz);
for n=1:length(invol);
    for m=1:6;
        temp=abs(fftshift(fft2(invol{n}(:,:,m),padsz,padsz)));
       % fftim = max(fftim,temp./max(temp(:)));
        fftim = max(fftim,temp./2^n);
    end;
end