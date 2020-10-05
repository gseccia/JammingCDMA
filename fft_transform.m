function [interval,frequency_sequence] = fft_transform(sequence, f)
    %Effettua la trasformata della sequenza e la restituisce insieme al
    %corretto intervallo di frequenze
    [nut,n_fft]=size(sequence);
    n_fft = n_fft*15;
    frequency_sequence=zeros(nut,n_fft);
    for i=1:nut
        frequency_sequence(i,:) = fftshift(abs(fft(sequence(i,:),n_fft)))/f;
    end
    interval = [-1/2:1/n_fft:1/2-1/n_fft]*f;
end