function [original_audio,v_t,fr_vec] = create_tx_signal(audio_paths)
    n = zeros(1, length(audio_paths));
        
    %Individua la lunghezza dei file audio presenti nei path specificati
    for i = 1:length(audio_paths)
        tmp = audioinfo(audio_paths{i});
        n(i) = tmp.TotalSamples;
    end
    
    %Definisce la lunghezza delle sequenze da considerare in base a quella
    %minore
    len = min(n)
    v_t = zeros(length(audio_paths), len*16);
    original_audio = zeros(length(audio_paths), len);
    fr_vec = zeros(1, length(audio_paths));
    
    %Legge i diversi file e ne ottiene una rappresentazione binaria
    for i = 1:length(audio_paths)
        [audio, audio_vec, fs] = from_audio_to_bits(audio_paths{i}, len);
        v_t(i, :) = audio_vec;
        fr_vec(i) = fs;
        original_audio(i,:) = audio;
    end
end

