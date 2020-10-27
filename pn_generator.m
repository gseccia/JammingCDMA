function [c_t] = pn_generator(Nuser,Lc,len_signal)
    %Genera una sequenza pseudo-casuale per ogni utente da utilizzare come
    %codice per la modulazione CDMA
    c_t_ref = zeros(Nuser,Lc);
    for i=1:Nuser
%         pn_seq_gen = comm.PNSequence('Polynomial', 'z^53 + z^17 + z^2 + 1','SamplesPerFrame',Lc,'InitialConditions',sequence_generator(1,53));
        pn_seq_gen = mseq(2,log2(Lc+1), Lc+1, i)';
        c_t_ref(i,:)=pn_seq_gen()';
    end
    
    %Replica tale codice per l'intera lunghezza del segnale
    c_t = repmat(c_t_ref,1,len_signal);
    
    %unipolar to bipolar conversion
    c_t(c_t==0)=-1;
end

