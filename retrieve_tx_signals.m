function [rec_sounds] = retrieve_tx_signals(rxbits, fr_vec, audio_paths)
    %Per ogni sequenza ricevuta ricostruisce e salva l'audio
    rec_sounds = zeros(size(rxbits, 1), size(rxbits, 2)/16);
    for i = 1:size(rxbits,1)
        file_name = strrep(audio_paths{i},'audio/','results/');
        [rx_audio] = from_bits_to_audio(file_name, fr_vec(i), rxbits(i,:));
        rec_sounds(i, :) = rx_audio;
    end
end

