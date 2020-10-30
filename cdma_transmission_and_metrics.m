% Implementazione di una comunicazione con protocollo d'accesso CDMA basato
% su Direct Sequence Spread Spectrum e modulazione BPSK, effettuata su 
% canale AWGN con eventuale presenza di jamming
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
    v_t_ext = sequence_extend(v_t_ref, Lc);
    v_t = sequence_extend(v_t_bpsk, Lc);
    
    % Spreading
    y_t = v_t.*c_t_ref;
    
    %Sovrapposizione dei segnali
    y_t_sum = sum(y_t,1);
        
    %In presenza di jamming
    if jammer_intensity_factor~=0
        %Calcolo della potenza del segnale trasmesso
        signal_power = rms(y_t(1,:))^2;
        %Calcolo della potenza del jammer partendo dalla potenza del
        %segnale
        Pj = signal_power*jammer_intensity_factor;

        %PI = Pj;
        t = [1:N]*Tc;
        
        % ------- Broadband Jammer ------------
        if jamming_type==0
%             t_picco = (0.5*rand()+0.25);
%             w = alpha*W;
%             jamming_signal = Pj*(sinc(w*(t-t_picco)));
            jammer = phased.BarrageJammer('ERP', 2*Pj, 'SamplesPerFrame', N);
            jamming_signal = real(jammer()');
%             jamming_signal =  wgn(1,N,Pj,'linear','real');
            
        % ------- Single Tone Jammer -------- 
        elseif jamming_type==1
            fb = (0.8*rand()+0.1)*W;
            %jamming_signal = sqrt(Pj)*ones(1, N);
            jamming_signal = sqrt(Pj)*cos(2*pi*fb*t);
                        
        % ------- Multi Tone Jammer --------     
        elseif jamming_type==2
            fb=[0.1:0.3:0.9]*W;
            jamming_signal = zeros(1, N);
            for i=1:length(fb)
                jamming_signal = jamming_signal + sqrt(Pj/length(fb))*cos(2*pi*fb(i)*t);
            end
        end
        
        Pj_calc = rms(jamming_signal)^2;
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
        
        Pj = 0;
        Pn = rms(y_awgn-y_t_sum)^2;
        PI = Pn;         
    end
    
    % Despreading
    Ps = rms(y_t(1,:))^2;
    v_des = y_awgn.*c_t_ref;
    Ps_des = rms(v_des(1,:))^2;
    
    % Integrazione
    v_des_integrated = integrate_signal(v_des,Lc);
    v_des_integrated_ext = sequence_extend(v_des_integrated, Lc);
        
    % Demodulazione BPSK
    rxbits = bpsk_demodulation(v_des_integrated);
    rxbits_ext = sequence_extend(rxbits, Lc);
    
    % Valutazione del BER per ogni utente
    [n_err, ber] = biterr(rxbits,v_t_ref,'row-wise');
    
    %Plot relativi alla comunicazione relativamente all'utente 1
    if plotting    
        %Potenze
        Potenza_jamming = Pj
        Potenza_rumore = Pn
        Potenza_interferenza = PI
        Potenza_segnale = Ps
        Potenza_despreading = Ps_des
        Potenza_interferenza_altri_utenti = Ps_des - PI - Ps
        % Dominio del tempo
        tc_interval = Tc:Tc:N*Tc;
        tb_interval = Tb:Tb:len_signal*Tb;
        % Dominio della frequenza        
        [tbf_interval, v_f] = fft_transform(v_t_ext, 1/Tc);
        [~, v_f_bpsk] = fft_transform(v_t, 1/Tc);
        [wf_interval, c_f] = fft_transform(c_t_ref, W);
        [tcf_interval, y_f] = fft_transform(y_t, 1/Tc);
        [~, y_f_sum] = fft_transform(y_t_sum, 1/Tc);
        [~, y_f_awgn] = fft_transform(y_awgn, 1/Tc);
        [~, v_des_f] = fft_transform(v_des, 1/Tc);
        [~, v_des_f_integrated] = fft_transform(v_des_integrated_ext, 1/Tc);
        [~, rxbits_f] = fft_transform(rxbits_ext, 1/Tc);
    
        
        figure
        %segnale in ingresso
        subplot(4,2,1)
        stairs(tc_interval, v_t_ext(1,:))
        ylim([-0.5,1.5])
        title("v(t)")

        subplot(4,2,2)
        plot(tbf_interval,v_f(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("v(f)")

        %segnale modulato BPSK
        subplot(4,2,3)
        stairs(tc_interval, v_t(1,:))
        ylim([-1.5,1.5])
        title("v_{BPSK}(t)")

        subplot(4,2,4)
        plot(tbf_interval, v_f_bpsk(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("v_{BPSK}(f)")

        %codice
        subplot(4,2,5)
        stairs(tc_interval,c_t_ref(1,:))
        ylim([-1.5,1.5])
        title("c(t)")

        subplot(4,2,6)
        plot(wf_interval,c_f(1,:))
        title("c(f)")

        %segnale spreaded
        subplot(4,2,7)
        stairs(tc_interval, y_t(1,:))
        ylim([-1.5,1.5])
        title("Spreading")

        subplot(4,2,8)
        plot(wf_interval, y_f(1,:))
        title("Spreading")

        figure
        %segnale complessivo
        subplot(3,2,1)
        stairs(tc_interval, y_t_sum)
        ylim([min(y_t_sum) * 1.2, max(y_t_sum) * 1.2])
        title("Segnale complessivo trasmesso")

        subplot(3,2,2)
        plot(tcf_interval, y_f_sum)
        title("Segnale complessivo trasmesso")


        % in presenza di jamming
        if jammer_intensity_factor~=0
            subplot(3,2,3)
            stairs(tc_interval, jamming_signal)
            ylim([min(jamming_signal) * 1.2 - 0.1, max(jamming_signal) * 1.2 + 0.1])
            title("Jamming Signal")

            [~, jamming_signal_f] = fft_transform(jamming_signal, 1/Tc);
            subplot(3,2,4)
            plot(tcf_interval,jamming_signal_f)
            title("Jamming Signal")
        end

        %segnale complessivo in ricezione
        subplot(3,2,5)
        stairs(tc_interval, y_awgn)
        ylim([min(y_awgn) * 1.2, max(y_awgn) * 1.2])
        title("Uscita del canale AWGN")

        subplot(3,2,6)
        plot(tcf_interval, y_f_awgn)
        title("Uscita del canale AWGN")


        figure
        %segnale despreaded
        subplot(3,2,1)
        stairs(tc_interval, v_des(1,:))
        ylim([min(v_des(1,:)) * 1.2, max(v_des(1,:)) * 1.2])
        title("Despreading")

        subplot(3,2,2)
        plot(tcf_interval, v_des_f(1,:))
        title("Despreading")

        %segnale integrato
        subplot(3,2,3)
        stairs(tb_interval, v_des_integrated(1,:))
        ylim([min(v_des_integrated(1,:)) * 1.2, max(v_des_integrated(1,:)) * 1.2])
        title("Integrazione")

        subplot(3,2,4)
        plot(tbf_interval, v_des_f_integrated(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("Integrazione")

        %demodulazione BPSK
        subplot(3,2,5)
        stairs(tb_interval, rxbits(1,:))
        ylim([-0.5,1.5])
        title("Demodulazione BPSK")

        subplot(3,2,6)
        plot(tbf_interval, rxbits_f(1,:))
        xlim([min(tcf_interval) max(tcf_interval)])
        title("Demodulazione BPSK")


        figure
        subplot(3,1,1)
        stairs(tb_interval, v_t_bpsk(1,:))
        hold on
        stairs(tb_interval, v_des_integrated(1,:))
        legend('BPSK','Integration')

        subplot(3,1,2)
        plot(tbf_interval, v_f_bpsk(1,:))
        hold on
        plot(tbf_interval, v_des_f_integrated(1,:))
        legend('BPSK','Integration')

        subplot(3,1,3)
        stairs(tb_interval, v_t_ref(1,:))
        hold on
        stairs(tb_interval, rxbits(1,:))
        ylim([-0.5,1.5])
        legend('Tx','Rx')
        
        % in presenza di jamming
        if jammer_intensity_factor~=0
            figure
            subplot(2,2,1)
            stairs(tc_interval, jamming_signal)
            ylim([-0.1 + min(jamming_signal) * 1.2, max(jamming_signal) * 1.2 + 0.1])
            title("Jamming Signal")

            [~, jamming_signal_f] = fft_transform(jamming_signal, 1/Tc);
            subplot(2,2,2)
            plot(tcf_interval,jamming_signal_f)
            title("Jamming Signal")
            
            jamming_signal_spreaded = jamming_signal.*c_t_ref(1,:);
            
            subplot(2,2,3)
            stairs(tc_interval, jamming_signal_spreaded)
            ylim([min(jamming_signal_spreaded) * 1.2, max(jamming_signal_spreaded) * 1.2])
            title("Jamming Signal Spreaded")

            [~, jamming_signal_f_spred] = fft_transform(jamming_signal_spreaded, 1/Tc);
            subplot(2,2,4)
            plot(tcf_interval,jamming_signal_f_spred)
            title("Jamming Signal")
        end    
    end
end