function [BER,SER] = LLL_ZF_SIC(n_symbols , n_iterations , bits_all , n_all , H_all , par)

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
        
        %Inicialização Erro
        err = 0;
        err_bits = 0;
        
        %parfor loop2=1:n_symbols
        for loop2 = 1:n_symbols
            % generate transmit symbol
            idx = bi2de(bits_all(:,:,loop2),'left-msb')+1;
            x = par.symbols(idx).'; 
            alpha = 2;
            s_real = [real(x);imag(x)];
            s_real = (s_real + ones(2*par.MT,1))/alpha;
            z_LR = pinv(T)*(s_real);
            
            % Sinal Recebido = Canal*SinalEnviado + Ruido
            n_real = [real(n_all(:,loop2,loop1));imag(n_all(:,loop2,loop1))];
            y_LR = H_LR*z_LR + sigma*n_real/alpha;
            
            %Lattice Estimação com ZF
            % Calculos necesarios para redução das lattice
            %% lattice reduced zero-forcing SIC (LR_ZF_SIC) detector
            [Q,R] = qr(H_LR);
            y_LR = Q'*y_LR;
            zhat(2*par.MT,1) = round( y_LR(2*par.MT,1)/R(2*par.MT,2*par.MT) ); 
            for ii = 2*par.MT-1:-1:1
                zhat(ii,1) =  round((y_LR(ii,1) - R(ii,ii+1:2*par.MT)*zhat(ii+1:2*par.MT,1)) / R(ii,ii));
            end
            xhat_quan = T*zhat;
            %Aproximar ao boundary
            r_temp = alpha*xhat_quan - ones(2*par.MT,1);
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
                for ii=1:L temp(ii)=(r_temp(kk)-V(ii))^2; end                
                indx=find(min(temp)==temp);
                r(kk,:)=V(indx(1));    
            end
            xhat_quan = r;
            xhat =  xhat_quan(1:par.MT,:) + 1i*xhat_quan((par.MT+1):2*par.MT,:);
            [a,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
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