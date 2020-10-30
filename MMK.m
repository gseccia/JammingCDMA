function [meanlen, len, leng] = MMK(lambda, mu,k,nv)

%arrive time stamp
at = zeros(1,nv);
%crossing time stamp
ct = zeros(1,nv); %tempo in cui inizia il servizio
%leave time stamp
lt = zeros(1,nv);

%%
gapt = exprnd(lambda,1,nv); %gap time (interarrivi)
st = exprnd(mu,1,nv); % service time del primo server
for i =2 : k
    st = [st; exprnd(mu,1,nv)]; % service time degli altri server
end

%%
%init the index of services
at(1) = 0;
for i =2 : nv
    at(i) = at(i-1) + gapt(i-1);
end

%the index number on the servers
index = zeros(1,k);
index(1) = 1;

%the pointer record the serviced number of each service station
pointer = zeros(1,k) + 1;

ct(1) = 0;
lt(1) = ct(1) + st(1,1);
pointer(1) = pointer(1) + 1;

%value and position of the packet which will leave
v = -1;
pos = -1;

%%
for i = 2 : k
    ct(i) = at(i);
    lt(i) = ct(i) + st(i,pointer(i));
    pointer(i) = pointer(i) + 1;
    index(i) = i;
end

%%
for i = (k+1):nv %per i restanti pacchetti
    
    
    [v, pos] = min(lt(index(1:k))); %prendi il primo server che si libera e l'istante in cui si libera
    
    if at(i) >= v % se il pacchetto arriva e il server è libero
        ct(i)  = at(i); %inizio del servizio all'arrival time
    else
        prev = index(pos); % mantengo l'informazione sulla prima macchina che si libera
        ct(i) = lt(prev); %il cross time sarà quando la macchina si libererà
    end
    
    lt(i) = ct(i) + st(pos,pointer(pos)); %registro il leave time
    pointer(pos) = pointer(pos)+1; %passo al prossimo service time
    index(pos) = i; %registro il pacchetto come ultimo servito per il server
end


%%
%legnth of the queue
lenat = zeros(nv,2); %per ogni arrival time qual è la lunghezza della coda in quell'istante
lenlt = zeros(nv,2); %per ogni leave time qual è la lunghezza della coda in quell'istante

for i = 1:nv %per ogni pacchetto
    lenat(i,1) = at(1,i); %inizializzo l'indice come il tempo di arrivo 
    lenat(i,2) = 1;  %inizializzo la lunghezza della coda come 1 (il primo)
    
    lenlt(i,1) = lt(1,i); %inizializzo l'indice come il tempo di partenza 
    lenlt(i,2) = -1; %inizializzo la lunghezza della coda come -1 (1 in servizio)
    
end
%%
len = [lenat; lenlt]; % associo le due code
len = sortrows(len,1); % ordino i tempi in ordine crescente


%%
leng = zeros(2,size(len,1)); %creo il vettore delle lunghezze di coda
available=k; %numero di server occupati
for i = 1:size(len,1) %all'inizio la coda è zero, a partire dal secondo
    leng(2,i) = len(i,1); 
    if len(i,2) > 0 % se è un arrivo
        if available > 0
            leng(1,i) = 0; %ho direttamente il servizio
            available = available -1;
        else
            leng(1,i) = leng(1,i-1) + 1; %aumento la coda
        end
    else %se è una partenza
        if leng(1,i-1) > 0  % e avevo elementi in coda
            leng(1,i) = leng(1,i-1)-1 ; %riduco la coda
            
        else %altrimenti
            leng(1,i) = 0; %ho solo consumato quelli in servizio
            available = available + 1;
        end %fine
    end %fine
    
    if available>k
        ME = MException('VerifyOutput:OutOfBounds', ...
             'Results are outside the allowable server limits');
        throw(ME);
    end
    if available<0
        ME = MException('VerifyOutput:OutOfBounds', ...
             'Results are outside the allowable server limits');
        throw(ME);
    end

end
%leng
%meanlen = mean(leng(1,:)); %fai la media delle lunghezze di coda

time1 = leng(2,:);
time2 = [0 leng(2,1:end-1)];
intervals = time1-time2;
meanlen = sum(leng(1,:).*intervals)/sum(intervals);


% %Vettore degli istanti temporali
% del = 0.1*1/lambda;
% Tmax = lt(end);
% tvec = [1:Tmax/del]*del;
% new_leng = zeros(1,length(tvec));
% % 
% leng_idx_prev = 1;
% for ii=2:Tmax/del
%     tc = (ii-1)*del; %istante temporale
%     while leng(2,leng_idx_prev) < tc
%         new_leng(ii-1) = leng(1,leng_idx_prev);
%         leng_idx_prev = leng_idx_prev+1;
%     end
% end
% meanlen = mean(new_leng);


end %fine


