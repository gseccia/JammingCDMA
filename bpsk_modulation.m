function [v_t] = bpsk_modulation(v_t_ref)
    %Modulazione BPSK utente per utente
    BPSKMod = comm.BPSKModulator();
    [Nuser,len_signal]=size(v_t_ref);
    v_t=zeros(Nuser,len_signal);
    
    for i=1:Nuser
        sing_sig_mod = real(BPSKMod(v_t_ref(i,:)'))';
        v_t(i,:)= sing_sig_mod;
    end
end

