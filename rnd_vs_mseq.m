c_mseq = pn_generator(1,511,10000,0);
c_pn = pn_generator(1,511,10000,1);

figure
subplot(2,1,1)
plot(autocorr(c_pn,'NumLags',3000))
subplot(2,1,2)
plot(autocorr(c_mseq,'NumLags',3000))

figure
plot(fftshift(abs(fft(c_pn,length(c_pn)*10))/length(c_pn)))
hold on
plot(fftshift(abs(fft(c_mseq,length(c_mseq)*10))/length(c_mseq)))
