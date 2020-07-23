function [BER,SER] = MODProcedure(n_symbols , n_iterations , bits_all , n_all , H_all , par)

%% PARAMETROS INICIAIS
%var_ruido = var(sqrt(N0)*n_real);  %Varianza Ruido
max_depth = 10; dp = 1;
theta1 = 0.8 ; theta2 = 4 ; theta3 = theta1*theta2 - 0.2;

%% COMEÇO SIMULAÇÃO
n_errors=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Simbolos
n_errors_bits=zeros(length(par.SNRdB_list),1); %Numero de Erros para cada SNR -> Bits

for loop1=1:n_iterations    
    Hc = H_all(:,:,loop1); %Matriz Canal com Distribuição Normal com media zero e varianza unitaria
    R = zeros(2*par.MT,1);  %Conjunto coordenadas confaveis
    H = [real(Hc) -imag(Hc) ; imag(Hc) real(Hc)];  

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
            y = H*s_real + sigma*n_real;
            
            %% Detecção
            %% Mod.1_Procedure
            %Input:{ y=Hx+n->reais ; H->reais ; max_depth -> max. iterações ; 
            % dp->depht recussion <= max_depth ; thetha1 ; thetha2 ; thetha3 }
        
            %Thetas são limitantes confiabilidade.
            %%% theta1 -> |zi' - ri'|<theta1 ;  0<theta1<1 -> smaller theta1, stronger
                % reliable condition.
            %%% theta2 -> var_noise*gii<theta2; Large values theta2 means coordinate with
                % large noise variance is realiable, theta2->covarianza ruido.
            %%% theta3 -> (var_noise*gii)*(|zi' - ri'|)<theta3; Melhora a condição de
                % confiabilidade é sempre menor que theta1*theta2.
            while 1
                T = LLL(H,par.MT); 
                H_tiu = H*T ; 
                %G = pinv(H_tiu'*H_tiu);
                r = pinv(H_tiu)*y;  % r = T^(-1)x + n = z + n ; onde z=T^(-1)x
                %Quantização Z
                z = round(r);
                alpha=2;
                f = z;
                for i=1:2*par.MT
                    if abs( z(i) - r(i) ) < theta1  %&&  var_ruido*G(i,i) < theta2  &&  (var_ruido*G(i,i))*(abs(z(i) - r(i))) < theta3
                        R(i,1) = i; %É Confiavel esta amostra
                    else
                        f(i,1) = 0; %É amostra não confiavel
                    end
                end %fim for        
                ylayer = y - H_tiu*f;
                for ii=1:size(R,1)
                    if R(ii,1) ~= 0; %O elemento de R não é NaN
                       Hlayer(:,ii) = H_tiu(:,ii);
                    end
                end
                if sum(abs(R))==0 || dp>max_depth
                    z_hat = pinv(Hlayer)*ylayer;
                    z_hat = round(1/alpha*z_hat-1/2*pinv(T)*ones(2*par.MT,1)); %Quantização
                    break;
                else %Continua o Algoritmo
                    y = ylayer; %Update para iteração
                    H = Hlayer; %Update para iteração
                    dp = dp + 1;
                end                
            end %fim WHILE
            xhat_quan = T*z_hat;
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
    


    
    
    
    
    
    
    
    

    
