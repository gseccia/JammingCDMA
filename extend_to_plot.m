function [interval,plot_sequence] = extend_to_plot(sequence, T)
    %Restituisce intervallo e sequenza utile per il plotting del segnale
    [nut,nlen]=size(sequence);
    plot_sequence=zeros(nut,nlen+1);
    for i=1:nut
        plot_sequence(i,:) = [sequence(i,:),sequence(i,end)];
    end
    interval = 0:T:T*nlen;
end