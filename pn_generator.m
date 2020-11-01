function [c_t] = pn_generator(Nuser,Lc,len_signal)
    %Genera una sequenza pseudo-casuale per ogni utente da utilizzare come
    %codice per la modulazione CDMA
    c_t_ref = zeros(Nuser,Lc);
    degree = log2(Lc+1);
    polys = primpoly(degree, 'all', 'nodisplay');
    k=0;
    for i=1:Nuser
        k = k+1;
        if k == length(polys)+1
            k = 1;
        end
        pn_seq_gen = comm.PNSequence('Polynomial', de2bi(polys(k)),'SamplesPerFrame',Lc,'InitialConditions',[zeros(1, degree-1) 1]);
        c_t_ref(i,:)=pn_seq_gen()';
    end
    
    %Replica tale codice per l'intera lunghezza del segnale
    c_t = repmat(c_t_ref,1,len_signal);
    
    %unipolar to bipolar conversion
    c_t(c_t==0)=-1;
end

