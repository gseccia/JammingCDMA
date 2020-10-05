function [rx_audio] = from_bits_to_audio(audio_name, fs, rxbits)
    %Legge una sequenza binaria e la trasforma in int16 prima di salvare su
    %un file audio 
    audio_uint_rec = zeros(1,length(rxbits)/16);
    j = 1;
    for i = 1:16:length(rxbits)
        audio_uint_rec(j) = bin2dec(num2str(rxbits(i:i+15)));
        j = j+1;
    end
    audio_uint_rec = uint16(audio_uint_rec);
    rx_audio = typecast(audio_uint_rec, 'int16');
    audiowrite(audio_name,rx_audio',fs);
end

