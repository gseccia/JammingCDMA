function [sequences] = sequence_generator(Nuser,length)
    %Definiti il numero di utenti e la lunghezza delle sequenze genera una
    %matrice di valori casuali dove ogni riga corrisponde a una sequenza
    %di bit lunga length relativa ad un particolare utente sul canale
    rng('shuffle')
    sequences=zeros(Nuser,length);
    for i=1:Nuser
        sequences(i,:) = randn(1,length)>0.0;
    end
end
