function [rxbits] = bpsk_demodulation(v_des_integrated)
    %Demodulazione BPSK utente per utente
    BPSKDemod = comm.BPSKDemodulator();
    [Nuser,len_signal] = size(v_des_integrated);
    rxbits=zeros(Nuser,len_signal);
    
    for i=1:Nuser
        sing_sig_demod = BPSKDemod(v_des_integrated(i,:)')';
        rxbits(i,:)= sing_sig_demod;
    end
end

