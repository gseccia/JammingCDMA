function [extended_sequences] = sequence_extend(sequence, Lc)
    %Ricevuta una matrice di sequenze in input estende ogni valore Lc volte
    [Nuser,length] = size(sequence);
    extended_sequences = zeros(Nuser,length*Lc);
    for i=1:Nuser
        extended_sequences(i,:) = rectpulse(sequence(i,:),Lc);
    end
end