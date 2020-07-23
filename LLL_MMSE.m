function [BER,SER] = LLL_MMSE(n_symbols , n_iterations , bits_all , n_all , H_all , par)

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations    
    Hc = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    H_real = [real(Hc) -imag(Hc) ; imag(Hc) real(Hc)];  
    T = LLL(H_real,par.MT); 
    H_LR = H_real*T;

    
    %parfor ind_db=1:length(par.SNRdB_list)        
    for ind_db = 1:length(par.SNRdB_list)
        sigma = sqrt(par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10)); %Sigma do Ruido AWGN
        %N0 = par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10);         
        Wh = pinv(H_LR'*H_LR + (sigma^2/par.Es)*(T')*T)*(H_LR');
        
        %Inicialização Erro
        err = 0;
        err_bits = 0;
        
        %parfor loop2=1:n_symbols
        for loop2 = 1:n_symbols
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).'; 
            s_real = [real(x);imag(x)];
            s_LR = pinv(T)*s_real;
            x_LR = H_LR*s_LR;

            % Sinal Recebido = Canal*SinalEnviado + Ruido
            n_real = [real(n_all(:,loop2,loop1));imag(n_all(:,loop2,loop1))];
            y_LR = x_LR + sigma*n_real;
            
            %Detecção Lattice com MMSE            
            xhat_real = Wh*y_LR; % inv(H)*y                        
            xhat_quan = MyQuan_LR(xhat_real,2*par.MT,T,length(par.symbols)); 
            xhat =  xhat_quan(1:par.MT,:) + 1i*xhat_quan((par.MT+1):2*par.MT,:);
            [none,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            bits_hat = par.bits(idxhat,:);
            %Conteo Erros
            err = err + sum(x~=xhat); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end %fim Loop2-> n_symbols
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end %fim ind_db->SNR
end%fim loop1

SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         