%clc; clear all
%par.MR = 4; % receive antennas 
%par.MT = 4; % transmit antennas (set not larger than MR!) 
%H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));

function T = LR_ELR(H,par)

%Inicio
C_tiu = ((H'*H)^(-1));
T_linha = eye(par.MR);


%% Algoritmo 1
while 1
    lambda = 0; %Reinicia lambda
    for ii=1:size(C_tiu,1) %Filas de C_tiu
        for kk=1:size(C_tiu,2) %Colunas de C_tiu 
            if ii==kk
               lambda(ii,kk) = 0;
            else
               lambda(ii,kk) = -round( C_tiu(ii,kk)/C_tiu(ii,ii) ); 
            end       
        end
    end
    Ckk = diag(C_tiu); %S� preceisso valores na diagonal
    [Largest_Ckk , k] = max(Ckk); %Pego a posi��o do max, como est� na diagonal, k-> [fila, coluna] do max. Largest   
    Delta = 0; %reiniciar Delta
    for i_tiu = 1:size(C_tiu,1) %PseudoCode Diz at� N., i_tiu diferente k
        if i_tiu ~= k
            del = -abs(lambda(i_tiu,k))^2*C_tiu(i_tiu,i_tiu) - conj(lambda(i_tiu,k))*C_tiu(i_tiu,k) - lambda(i_tiu,k)*C_tiu(k,i_tiu);
            if del>=0  %Revisar esta condi��o
                Delta(i_tiu) = -abs(lambda(i_tiu,k))^2*C_tiu(i_tiu,i_tiu) - conj(lambda(i_tiu,k))*C_tiu(i_tiu,k) - lambda(i_tiu,k)*C_tiu(k,i_tiu);
            else
                Delta(i_tiu) = NaN;
            end
        end
    end    
    if sum(Delta) == 0
       break; 
    end    
    [Delta_max , i] = max(Delta);    
    %Update
    T_linha(:,k) = T_linha(:,k) + lambda(i,k)*T_linha(:,i);
    C_tiu(:,k) = C_tiu(:,k) + lambda(i,k)*C_tiu(:,i); %Coluna C_tiu
    C_tiu(k,:) = C_tiu(k,:) + conj(lambda(i,k))*C_tiu(i,:); %Fila C_tiu
end %FIM WHILE

T = (T_linha^(-1))';
%H_tiu = H*T;


