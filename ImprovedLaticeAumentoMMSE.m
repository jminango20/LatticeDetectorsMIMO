%% Improved Lattice Reduction Aided Detections MMSE for MIMO Systems
function xhat_quan = ImprovedLaticeAumentoMMSE(Hreal,par,xreal,w)  %Hc e xc ->Canal e Simbolo Complexo 
    sigma = sqrt(var(w));  %Desvio Padrão
    alpha = 2; %Alfabeto Lattice
    x = (1/alpha)*xreal - (1/2)*ones(2*par.MT,1);  %Transformo a Lattice pertenece a Z
    H = alpha*Hreal;
    Hex = [H;sigma*eye(2*par.MR)];
    
    %Sinal Transmitido em Dimensões de Lattices
    y = H*x + w; %1 Jeito  
    yex = [y;zeros(2*par.MR,1)];
    
    %1Paso
    Tex = LLL(Hex,size(Hex,2)/2);  %Uso LLL para Hex->MatrizExtendida
    H_tiu_ex = Hex*Tex;
    [Q_tiu_ex,R_tiu_ex]=qr(H_tiu_ex,0);
    e = abs(R_tiu_ex(2*par.MT,2*par.MT));
    %2Paso Aumento Matriz
    Ha = [Hex yex; zeros(1,2*par.MT) e];
    Qa = [Q_tiu_ex zeros(2*par.MT+2*par.MR,1) ; zeros(1,2*par.MT) 1];
    Ra = [R_tiu_ex Q_tiu_ex'*yex ; zeros(1,2*par.MT) e];
    Ta = [Tex zeros(2*par.MT,1) ; zeros(1,2*par.MT) 1];
    
    
    %% LLL Aumentado
    Q_tiu_a = Qa; R_tiu_a = Ra; T_tiu_a = Ta; m = 2*par.MT + 1; 
    for k = m-1:-1:1
        u = round( R_tiu_a(k,m) / R_tiu_a(k,k) );
        R_tiu_a(1:k,m) = R_tiu_a(1:k,m) - u*R_tiu_a(1:k,k);
        T_tiu_a(:,m) = T_tiu_a(:,m) - u*T_tiu_a(:,k);
    end
    H_tiu_a = Ha*T_tiu_a;
    
    %Obtenho os pontos da lattice
    xhat_LR = (-1)*T_tiu_a( 1:2*par.MT , 2*par.MT+1 );               
    xrealhat = (alpha)*xhat_LR + (alpha/2)*ones(2*par.MT,1);  %Transformo a Lattice pertenece a Z
    
    %Quantização ao Alfabeto
    r_temp = xrealhat;
    switch length(par.symbols)
         case 2
             V=[1;-1];
         case 4
             V=[1;-1];        
         case 8
             V=[1;-1;3;-3];
         case 16
             V=[1;-1;3;-3];
         case 32
             V=[1;-1;3;-3;5;-5];       
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
    xhat_quan = r;