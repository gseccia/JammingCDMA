function [v_t_int] = integrate_signal(v_t_des,Lc)
    %Ricevuto il segnale despreaded si integra l'informazione trasmessa
    [Nuser,length]=size(v_t_des);
    len_signal = length/Lc;     %Lunghezza del segnale originale trasmesso
    v_t_int=zeros(Nuser,len_signal);

    %Integrale tra n*Tb e (n+1)*Tb del segnale despreaded per ogni utente
    for k=1:Nuser
        for i=1:len_signal
            for j=(i-1)*Lc+1:i*Lc
                v_t_int(k,i)=v_t_int(k,i)+v_t_des(k,j);
            end
        end
    end
    v_t_int = v_t_int./Lc;
end