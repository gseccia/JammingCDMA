function [original_audio, audio_vec, fs] = from_audio_to_bits(audio_path, len)
    %Legge un file audio, lo interpreta come uint16 e lo converte in una
    %sequenza binaria
    [original_audio, fs]=audioread(audio_path, 'native');
    original_audio = original_audio(1:len)';
    audio = original_audio;
    audio_uint = typecast(audio, 'uint16');

    audio_vec = zeros(1,length(audio_uint)*16);
    j=1;
    for i=1:length(audio_uint)
        audio_vec(j:j+15) = dec2bin(audio_uint(i),16) - '0';
        j=j+16;
    end
end