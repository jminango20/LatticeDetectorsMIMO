function xhat_quan = sampling(H_real,par,xreal,w)
    T = LLL(H_real,par.MT);
    H_tiu = H_real*T;
    [Q,R] = qr(H_tiu);
    %Para a qr do Matlab
    dR = 2*(diag(R) >= 0) - 1; % avoid using sign which can return 0, destroying qr-Diag(R) é positiva
    nqr = length(dR);
    Q(:,1:nqr) = bsxfun(@times, Q(:,1:nqr), dR'); R(1:nqr,:) = bsxfun(@times,R(1:nqr,:), dR);

    %Começa o Algoritmo
    n = size(H_tiu,2); %Já é 2n para o calculo de Kcal    
    alpha=2;
      
    yreal = H_real*xreal + w;  %USAR AO FINAL
    
    %% Começa a Simulação    
    x = (xreal + ones(2*par.MT,1))/alpha;  %Transformo a Lattice pertenece a Z
    z = pinv(T)*x;
    y = H_tiu*z + w/alpha;
    %Sinal Transmitido em Dimensões de Lattices
    y = Q'*y;
    
    %Achar rho em função do K
    K = 2; %É fixado
    rho = 130;   %rho>n então K<(exp*n)^2  ;  rho=n então K=(exp*n)^2 
%     while(1)        
%         Kcal = (exp(1)*rho)^(n/rho);
%         if Kcal>=(K-0.5) && Kcal<=(K+0.5)
%             break;
%         end
%         rho = rho + 1;
%     end    
    
    [minimo,pos] = min(diag(R));
    r_x = 0; z_tiu = zeros(n,1);
    
    for k=1:K %Iterações
    
        for ii=size(y,1):-1:1 %Rand_SIC_A 
            %c(ii) = ((log(rho))/(minimo))*(R(ii,ii)^2);
            c(ii) = 2;
            for jj = ii+1:size(y,1)
                r_x(jj) = R(ii,jj)*z_tiu(jj,1);  
            end
            Interferencia = sum(r_x);
            z_tiu(ii,1) = ( y(ii,1) - Interferencia ) / (R(ii,ii));            
            zhat = PDF(z_tiu(ii,1),c(ii));  
            z_tiu(ii,1) = zhat;    
        end %Fim ii        
    zhatTodos(:,k)=z_tiu;  
    z_tiu = 0; Interferencia = 0; r_x=0;
    end %Fim Iterações
    

    for ii=1:size(zhatTodos,2)
        r_temp = 2*(T*zhatTodos(:,ii)) - ones(2*par.MT,1);
        
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
            for iii=1:L
                temp(iii)=(r_temp(kk)-V(iii))^2;                                
            end            
            indx=find(min(temp)==temp);            
            r(kk,:)=V(indx(1));            
        end        
        
        xhat(:,ii) = r;
    end %fim For
    
    yoriginal = repmat(yreal,1,size(xhat,2));
    [minimo,pos] = min(sum(abs(yoriginal - H_real*xhat).^2));
    xhat_quan = xhat(:,pos);        
    

   