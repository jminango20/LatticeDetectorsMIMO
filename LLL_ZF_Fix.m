function [BER,SER] = LLL_ZF_Fix(n_symbols , n_iterations , bits_all , n_all , H_all , par)

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations    
    Hc = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    T_fixed = MLLL(Hc,par.MT);
    H_fixed = Hc*T_fixed;

    
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
            s_fixed = pinv(T_fixed)*x;
            x_fixed = H_fixed*s_fixed;
            % Sinal Recebido = Canal*SinalEnviado + Ruido
            y_LR = x_fixed + sigma*n_all(:,loop2,loop1);

            %Lattice Estimação com ZF-fix Complexo
            % Calculos necesarios para redução das lattice
            %% lattice reduced zero-forcing FIX (LR_ZF_FIX) detector
            xhat_col = H_fixed\y_LR; % inv(H)*y
            %Quantization for lattice reduction, quantization is performed according
            %to the constellation points of inv(T)*s
            xhat_quan = MyQuan_CLR(xhat_col,2*par.MT,T_fixed,length(par.symbols)); 
            %Converting from real to complex
            xhat =  xhat_quan(1:par.MT,:) + 1i*xhat_quan((par.MT+1):2*par.MT,:);
            %Quantization and generating the index of the bits
            [a,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            %Finding the corresponding bits according to the indexing
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