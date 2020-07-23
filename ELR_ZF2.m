function [BER,SER] = ELR_ZF2(n_symbols , n_iterations , bits_all , n_all , H_all , par) %ELR Element lattice reduced zero-forcing detector.

n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations    
    H = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    
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
            xtx = par.symbols(idx).'; 
            %NOVO 
            alpha = 2;
            x = (xtx + ( ones(par.MT,1) + 1i*ones(par.MT,1) ) )/alpha;  
            % Calculos necesarios para redução das lattice
            T_ELR = LR_ELR_2(H,par);
            z_ELR = pinv(T_ELR)*x;
            H_ELR = H*T_ELR;  
            y = H_ELR*z_ELR + (sigma*n_all(:,loop2,loop1))/alpha;  

            %% ELR Element lattice reduced zero-forcing detector.  
            y = pinv(H_ELR)*y;
            zhat = round(y);
            shat = T_ELR*zhat;
            shat = alpha*shat - ( ones(par.MT,1) + 1i*ones(par.MT,1) );   
            r_temp = [real(shat);imag(shat)];
            %Quantization and generating the index of the bits
            switch length(par.symbols)
                 case 2
                      V=[1;-1];
                 case 4
                      V=[1;-1];
                 case 16
                      V=[1;-1;3;-3];
                 case 64
                      V=[1;-1;3;-3;5;-5;7;-7]; 
            end
            L=length(V);   
            for kk=1:length(r_temp)            
                for ii=1:L
                    temp(ii)=(r_temp(kk)-V(ii))^2;                          
                end                
                indx=find(min(temp)==temp);
                r(kk,:)=V(indx(1));      
            end
            s_hat = r(1:par.MT,1)+1i*r(par.MT+1:end,1);  
            [a,idxhat] = min(abs(s_hat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
            %Finding the corresponding bits according to the indexing
            xhat = par.symbols(idxhat).'; %Decodificação do x_mmse - Quantização
            bits_hat = par.bits(idxhat,:);            
            %Conteo Erros
            err = err + sum(xtx~=xhat); %Erro De Simbolos
            err_bits = err_bits + sum(sum(bits_all(:,:,loop2)~=bits_hat)); %Erro De Bits
        end %fim Loop2-> n_symbols
        n_errors(ind_db)=err;
        n_errors_bits(ind_db)=err_bits;
    end %fim ind_db->SNR
end%fim loop1

SER = n_errors/(par.MT*n_iterations*n_symbols); %SER         
BER = n_errors_bits/(par.MT*par.Q*n_iterations*n_symbols); %BER         