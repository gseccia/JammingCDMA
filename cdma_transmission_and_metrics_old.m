
function [n_err, ber, rxbits] = cdma_transmission_and_metrics(v_t_ref,Lc,Tb,EbNo, plotting,jammer_intensity_factor,jamming_type,alpha)
    % Matrice dei segnali, una riga per ogni utente
    [Nuser,len_signal]=size(v_t_ref);
    
    Tc = Tb/Lc; % Chip interval
    N=len_signal*Lc; % Lunghezza della sequenza codificata
    
    EbNo = db(EbNo, 'power');        
    W = 1/Tc;
    
    % Modulazione BPSK
    v_t_bpsk = bpsk_modulation(v_t_ref);

    % Generazione codici PN
    c_t_ref = pn_generator(Nuser,Lc,len_signal);
    
    % Estensione dei segnali per renderli confrontabili ai codici
    v_t = sequence_extend(v_t_bpsk, Lc);
    
    % Spreading
    y_t = v_t.*c_t_ref;
      
    %Sovrapposizione dei segnali
    y_t_sum = sum(y_t,1);
        
    %In presenza di jamming
    if jammer_intensity_factor~=0
        %Calcolo della potenza del segnale trasmesso
        signal_power = rms(y_t_sum)^2;
        %Calcolo della potenza del jammer partendo dalla potenza del
        %segnale
        Pj = signal_power*jammer_intensity_factor
        t = [1:N]*Tc;
        
        % ------- BROADBAND JAMMER ------------
        if jamming_type==0
            t_iniziale = 0.25*N*Tc;
            w = alpha*W;
            
            jamming_signal = (Pj)*(sinc(w*(t-t_iniziale)));
          
        % ------- Single Tone Jammer -------- 
        elseif jamming_type==1
            fb = (0.4*rand()+0.1)*W; % scelta casualmente
            jamming_signal = sqrt(2*Pj)*cos(2*pi*fb*t);
                        
        % ------- Multi Tone Jammer --------     
        elseif jamming_type==2
            fb=0.1*[1:10]*W;
            jamming_signal = zeros(1,N);
            for i=1:length(fb)
                jamming_signal = jamming_signal + sqrt(2*Pj/length(fb))*cos(2*pi*fb(i)*t);
            end
            
        end
        %Definizione del canale AWGN con Eb/N0 fissato rispetto alla
        %potenza del segnale da trasmettere (informazione+jammer)
        channel = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)','EbNo',EbNo,'SignalPower',rms(y_t_sum+jamming_signal)^2);
        y_awgn = channel(jamming_signal+y_t_sum);
        
        Pn = rms(y_awgn-jamming_signal-y_t_sum)^2;
        PI = Pj + Pn;
    else %In assenza di jamming
        %Definizione del canale AWGN con Eb/N0 fissato rispetto alla
        %potenza del segnale da trasmettere
        channel = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)','EbNo',EbNo,'SignalPower',rms(y_t_sum)^2);
        y_awgn = channel(y_t_sum);
        
        Pn = rms(y_awgn-y_t_sum)^2;
        PI = Pn;
    end
    
    % Despreading
    Ps = rms(y_t(1,:))^2;
    v_des = y_awgn.*c_t_ref;
    Ps_des = rms(v_des(1,:))^2;
    
    % Integrazione
    v_des_integrated = integrate_signal(v_des,Lc);   
        
    % Demodulazione BPSK
    rxbits = bpsk_demodulation(v_des_integrated);
    
    % Valutazione del BER per ogni utente
    [n_err, ber] = biterr(rxbits,v_t_ref,'row-wise');
    
    %Plot relativi alla comunicazione relativamente all'utente 1
    if plotting
        % Dominio del tempo
        tc_interval = 0:Tc:N*Tc;
        [tb_interval, v_t_plot] = extend_to_plot(v_t_ref, Tb);
        [~, v_t_bpsk_plot] = extend_to_plot(v_t_bpsk, Tb);
        [~, c_t_plot] = extend_to_plot(c_t_ref, Tc);
        [~, y_t_plot] = extend_to_plot(y_t, Tc);
        [~, y_t_sum_plot] = extend_to_plot(y_t_sum, Tc);
        [~, y_awgn_plot] = extend_to_plot(y_awgn, Tc);
        [~, v_des_plot] = extend_to_plot(v_des, Tc);
        [~, v_des_integrated_plot] = extend_to_plot(v_des_integrated, Tb);
        [~, rxbits_plot] = extend_to_plot(rxbits, Tb);
        
        figure
        %segnale in ingresso
        subplot(4,3,1)
        stairs(tb_interval, v_t_plot(1,:))
        ylim([-0.5,1.5])
        title("v(t)")
        
        %segnale modulato BPSK
        subplot(4,3,4)
        stairs(tb_interval, v_t_bpsk_plot(1,:))
        ylim([-1.5,1.5])
        title("Modulazione BPSK")
        
        %codice
        subplot(4,3,7)
        stairs(tc_interval,c_t_plot(1,:))
        ylim([-1.5,1.5])
        title("c(t)")
        
        %segnale spreaded
        subplot(4,3,10)
        stairs(tc_interval, y_t_plot(1,:))
        ylim([-1.5,1.5])
        title("Spreading")

        %segnale complessivo
        subplot(4,3,5)
        stairs(tc_interval, y_t_sum_plot)
        ylim([min(y_t_sum_plot) * 1.2, max(y_t_sum_plot) * 1.2])
        title("Segnale complessivo trasmesso")
        
        %segnale complessivo in ricezione
        subplot(4,3,11)
        stairs(tc_interval, y_awgn_plot)
        ylim([min(y_awgn_plot) * 1.2, max(y_awgn_plot) * 1.2])
        title("Uscita del canale AWGN")
        
        %segnale despreaded
        subplot(4,3,3)
        stairs(tc_interval, v_des_plot(1,:))
        ylim([min(v_des_plot(1,:)) * 1.2, max(v_des_plot(1,:)) * 1.2])
        title("Despreading")
        
        %segnale integrato
        subplot(4,3,6)
        stairs(tb_interval, v_des_integrated_plot(1,:))
        ylim([min(v_des_integrated_plot(1,:)) * 1.2, max(v_des_integrated_plot(1,:)) * 1.2])
        title("Integrazione")
        
        %demodulazione BPSK
        subplot(4,3,9)
        stairs(tb_interval, rxbits_plot(1,:))
        ylim([-0.5,1.5])
        title("Demodulazione BPSK")
        
        %Bit Error Rate
        subplot(4,3,12)
        text(0.2,0.5,'BER: ')
        text(0.5,0.5,num2str(ber(1)))
        axis off
        
        % in presenza di jamming
        if jammer_intensity_factor~=0
            [~, jamming_signal_plot] = extend_to_plot(jamming_signal, Tc);
            subplot(4,3,8)
            stairs(tc_interval, jamming_signal_plot)
            title("Jamming Signal")
        end
        
        sgtitle('Analisi nel Dominio del Tempo per Utente 1')
        % Dominio della frequenza        
        [tbf_interval, v_f_plot] = fft_transform(v_t_ref, 1/Tb);
        [~, v_f_bpsk_plot] = fft_transform(v_t_bpsk, 1/Tb);
        [wf_interval, c_f_plot] = fft_transform(c_t_ref, W);
        [tcf_interval, y_f_plot] = fft_transform(y_t, 1/Tc);
        [~, y_f_sum_plot] = fft_transform(y_t_sum, 1/Tc);
        [~, y_f_awgn_plot] = fft_transform(y_awgn, 1/Tc);
        [~, v_des_f_plot] = fft_transform(v_des, 1/Tc);
        [~, v_des_f_integrated] = fft_transform(v_des_integrated, 1/Tb);
        [~, rxbits_f_plot] = fft_transform(rxbits, 1/Tb);
     
        figure
        
        %segnale in ingresso
        subplot(4,3,1)
        plot(tbf_interval,v_f_plot(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("v(f)")
                
        %segnale modulato BPSK
        subplot(4,3,4)
        plot(tbf_interval, v_f_bpsk_plot(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("Modulazione BPSK")
        
        %codice
        subplot(4,3,7)
        plot(wf_interval,c_f_plot(1,:))
        title("c(f)")

        %segnale spreaded
        subplot(4,3,10)
        plot(wf_interval, y_f_plot(1,:))
        title("Spreading")
        
        %segnale complessivo
        subplot(4,3,5)
        plot(tcf_interval, y_f_sum_plot)
        title("Segnale complessivo trasmesso")

        %segnale complessivo in ricezione
        subplot(4,3,11)
        plot(tcf_interval, y_f_awgn_plot)
        title("Uscita del canale AWGN")
      
        %segnale despreaded
        subplot(4,3,3)
        plot(tcf_interval, v_des_f_plot(1,:))
        title("Despreading")
        
        %segnale integrato
        subplot(4,3,6)
        plot(tbf_interval, v_des_f_integrated(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("Integrazione")
        
        %demodulazione BPSK
        subplot(4,3,9)
        plot(tbf_interval, rxbits_f_plot(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("Demodulazione BPSK")
        
        % in presenza di jamming
        if jammer_intensity_factor~=0
            [~, jamming_signal_f_plot] = fft_transform(jamming_signal, 1/Tc);
            subplot(4,3,8)
            plot(tcf_interval,jamming_signal_f_plot)
            title("Jamming Signal")
        end
        
        figure
        subplot(3,1,1)
        stairs(tb_interval, v_t_bpsk(1,:))
        hold on
        stairs(tb_interval, v_des_integrated(1,:))
        legend('BPSK','Integration')

        subplot(3,1,2)
        plot(tbf_interval, v_f_bpsk_plot(1,:))
        hold on
        plot(tbf_interval, v_des_f_integrated(1,:))
        legend('BPSK','Integration')

        subplot(3,1,3)
        stairs(tb_interval, v_t_ref(1,:))
        hold on
        stairs(tb_interval, rxbits(1,:))
        legend('Tx','Rx')
    end
end