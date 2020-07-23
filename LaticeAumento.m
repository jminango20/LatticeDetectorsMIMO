%% Lattice Aumentado
function xhat_quan = LaticeAumento(H_real,par,s_real,w)
    T = LLL(H_real,par.MT);
    H_tiu = H_real*T;
    [Q,R] = qr(H_tiu);
    m = size(H_real,2); n = size(H_real,1); %Coluna m=2Nt e Fila n=2Nr
    delta = 3/4; alfa = 1/(delta*1/4);    %delta->parametro LLL - se mudo, mudar em LLL também
    %Calculo Parametro Otimização t
    [aH_tiu, pos] = min(abs(diag(R)));  
    %e = 1/(2*sqrt(2)*alfa^((m-1)/2));  %e<=1/(2*sqrt(2)*alfa^((m-1)/2)) por definição
    e = 2^(-m/4);
    t = e*aH_tiu;
    %t = R(m,m)+0.5;  %Para o IMPROVED
    
    %Envio Sinais        
    y = H_real*s_real + w;
    %y = (1/2)*(y + H_real*ones(n,1)); %Para que pertenezca à lattice - Escalamento e Desplazamento
    
    %Aumento Matriz H
    Haum_tiu = [H_real -y ; zeros(1,m) t];
    Taum_tiu = LLL(Haum_tiu,size(Haum_tiu,2)/2);   %LLL muliplica por 2, por isso faço aqui a divissão.
    %Haumrea_tiu = Haum_tiu*Taum_tiu;             %Podemos leer rapidamente a menssagem de Taum_tiu=[+-x;+-1]-na primeira coluna de Taum_tiu    
    
    %% FORMA 1: Pego só a Primeira Coluna
%    z_hat = Taum_tiu(1:m,1)/Taum_tiu(m+1,1);
%    z_hat = round(z_hat);                        %Faço Quantização
%    xhat_quan = z_hat;           %Volto ao original
    
    %% FORMA 2: Para melhorar o desempenho uso todas as colunas de Taum_tiu
     posi = ((abs(Taum_tiu(m+1,:)))==1);
     if sum(posi) ~= 0
        k = size(Taum_tiu(:,posi),2);
     else
         k = 1;
         posi = 1;
     end
     Taum_pos = Taum_tiu(:,posi);
     for ii = 1:k
         zhatTodos(:,ii) = round(Taum_pos(1:m,ii)/Taum_pos(m+1,ii));
     end       
     ytodos = repmat(y,1,k);
     [minimo,posicion] = min(sum(abs(ytodos - H_real*zhatTodos).^2));
     xhat_quan = zhatTodos(:,posicion);        
     
    