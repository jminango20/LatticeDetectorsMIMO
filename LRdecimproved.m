function [BER,SER] = LRdecimproved(n_symbols , n_iterations , bits_all , n_all , H_all , par)

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations    
    Hc = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    H_real = [real(Hc) -imag(Hc) ; imag(Hc) real(Hc)];  
    
    %parfor ind_db=1:length(par.SNRdB_list)        
    for ind_db = 1:length(par.SNRdB_list)
        sigma = sqrt(par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10)); %Sigma do Ruido AWGN
        %N0 = par.MT*par.Es*10^(-par.SNRdB_list(ind_db)/10);         
        
        %Inicialização Erro
        err = 0;
        err_bits = 0;
        
        %parfor loop2=1:n_symbols
        for loop2 = 1:n_symbols
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).'; 
            s_real = [real(x);imag(x)];
            
            % Sinal Recebido = Canal*SinalEnviado + Ruido
            n_real = [real(n_all(:,loop2,loop1));imag(n_all(:,loop2,loop1))];
            w = sigma*n_real;

            %% JUAN - From LR Aided Toward MLD in MIMO - Só um Elemento            
            xhat_real = LRdec_improved(par,s_real,H_real,w);
            %xhat_quan = LRdec_improved_16QAM(y_real,H_real,par);
            %xhat_quan = LRdec_improved_16QAM_coset(y_real,H_real,par);            
            
            xhat =  xhat_real(1:par.MT,:) + 1i*xhat_real((par.MT+1):2*par.MT,:);
            [a,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            x_hat = par.symbols(idxhat).'; %Decodificação do x_zf - Quantização
            bits_hat = par.bits(idxhat,:);   

            %Conteo Erros
            err = err + sum(x~=x_hat); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end %fim Loop2-> n_symbols
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end %fim ind_db->SNR
end%fim loop1

SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         