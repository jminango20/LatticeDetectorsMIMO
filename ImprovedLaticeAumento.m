%% Improved Lattice Reduction Aided Detections for MIMO Systems
function xhat_quan = ImprovedLaticeAumento(Hreal,par,xreal,w,opcion)  %Hc e xc ->Canal e Simbolo Complexo 

   alpha=2;
   
   %% Geração da LATTICE de ACORDO ao ALFABETO USADO
   %decimal = 0:2^(par.MT*log2(par.M))-1;                   %Geração todos os possiveis número decimais 
   %binary = de2bi(decimal(:) ,'left-msb')';    %Conversão de Binario a Decimal
   %hmodem = modem.qammod('M',par.M, 'SymbolOrder', 'Gray', 'InputType', 'bit' );       %Modulador QAM
   %all_sym_mimo = modulate(hmodem,binary);        %Todas as combinações entre Mt e M (olho)
   %all_sym_mimo_real = [real(all_sym_mimo);imag(all_sym_mimo)];   
   %Lattice = (1/alpha)*all_sym_mimo_real - (1/2)*ones(2*par.MT,size(all_sym_mimo_real,2));
   
   
   %% Geração da LATTICE de ACORDO ao ALFABETO USADO - Isto é Suficente
   simbolos = [real(par.symbols);imag(par.symbols)];
   Lattice = (1/alpha)*simbolos - (1/2)*ones(size(simbolos));
   
    %% Começa a Simulação
    x = (1/alpha)*xreal - (1/2)*ones(2*par.MT,1);  %Transformo a Lattice pertenece a Z
    H = alpha*Hreal;
    
    %Sinal Transmitido em Dimensões de Lattices
    y = H*x + w; %1 Jeito    
    %yr = Hreal*xreal + w;
    %y =  yr - (1/2)*H*ones(2*par.MT,1) %2Jeito as duas linhas
    
    %Algoritmo LLL
    T = LLL(H,par.MT);  %Dá o mesmo fazer como H e com Hreal
    H_tiu = H*T;
   
    %ZF de H_tiu
    yhat = pinv(H_tiu)*y;
    zhat = T * round(yhat);  %round(yhat) -> Quantização em z
    xLR = zhat;
    %Quantização ao Alfabeto
    %xLR = alpha*T*zhat+alpha/2*ones(2*par.MT,1);
    
    %% DIFERENTES
        switch opcion
            case 1 %COMO DIZ O PAPER - Não Olho se Pertenecem à Lattice
                %Saber se xLR pertence à Lattice 
                i = ismember(xLR,Lattice);  %Se xLR ismeber Lattice i=1 / Se xLR notismeber Lattice i=0
                if (sum(i))~=size(xLR,1)  %Há elemento que não pertenecem à Lattice
                    for ii = 1:size(T,2)
                        xdelta(:,2*ii-1) = xLR + T(:,ii);
                        xdelta(:,2*ii) = xLR - T(:,ii);
                    end
                    yk = repmat(y,1,size(xdelta,2));
                    [a,posicion] = min(sum(abs(yk - H*xdelta).^2));
                    xhatLR = xdelta(:,posicion);
                else
                    xhatLR = xLR;
                end
                r_temp = alpha*xhatLR + alpha/2*ones(2*par.MT,1);
               %r_temp = alpha*T*zhat+alpha/2*ones(2*par.MT,1);
               
            case 2
                %Saber se xLR pertence à Lattice  Olho Primeiramente se Pertenecem à Lattice
                i = ismember(xLR,Lattice);  %Se xLR ismeber Lattice i=1 / Se xLR notismeber Lattice i=0
                if (sum(i))~=size(xLR,1)  %Há elemento que não pertenecem à Lattice
                    [pos,valor] = find(i==0);  %Vou encontrar quais simbolos do xLR estão fora da Lattice
                    [pos1,valor] = find(i==1); %Vou encontrar quais simbolos do xLR estão na Lattice
                    for ii = 1:size(T,2)
                        xdelta( pos , 2*ii-1 ) = xLR(pos) + T(pos,ii);
                        xdelta( pos1 , 2*ii-1 ) = xLR(pos1);
                        xdelta( pos , 2*ii ) = xLR(pos) - T(pos,ii);
                        xdelta( pos1 , 2*ii ) = xLR(pos1);                       
                    end
                    yk = repmat(y,1,size(xdelta,2));
                    [a,posicion] = min(sum(abs(yk - H*xdelta).^2));
                    xhatLR = xdelta(:,posicion);
                else
                    xhatLR = xLR;
                end
               r_temp = alpha*xhatLR + alpha/2*ones(2*par.MT,1);
               %r_temp = alpha*T*zhat+alpha/2*ones(2*par.MT,1);

          
            case 3
               %ZF de H_tiu
               theta1 = 0.2;
               
               yhat = pinv(H_tiu)*y;
               zhat = round(yhat);  %round(yhat) -> Quantização em z

               zconf = NaN(size(zhat,1),1);        %zconfiavel inizializo
               Tconf = NaN(size(T));               %Tconfiavel inizializo
               Hconf = NaN(size(H_tiu));           %Hconfiavel inizializo

               
               [confiavel , valor] = find(abs( yhat - zhat )<theta1);
               
                   
               %Quantização ao Alfabeto
               %xLR = alpha*T*zhat+alpha/2*ones(2*par.MT,1);

                %Saber se xLR pertence à Lattice 
                i = ismember(xLR,Lattice);  %Se xLR ismeber Lattice i=1 / Se xLR notismeber Lattice i=0
                if (sum(i))~=size(xLR,1)  %Há elemento que não pertenecem à Lattice
                    [pos,valor] = find(i==0);  %Vou encontrar quais simbolos do xLR estão fora da Lattice
                    [pos1,valor] = find(i==1); %Vou encontrar quais simbolos do xLR estão na Lattice
                    for ii = 1:size(T,2)
                        xdelta( pos , 2*ii-1 ) = xLR(pos) + T(pos,ii);
                        xdelta( pos1 , 2*ii-1 ) = xLR(pos1);
                        xdelta( pos , 2*ii ) = xLR(pos) - T(pos,ii);
                        xdelta( pos1 , 2*ii ) = xLR(pos1);                       
                    end
                    yk = repmat(y,1,size(xdelta,2));
                    [a,posicion] = min(sum(abs(yk - H*xdelta).^2));
                    xhatLR = xdelta(:,posicion);
                else
                    xhatLR = xLR;
                end
               r_temp = alpha*xhatLR + alpha/2*ones(2*par.MT,1);
               %r_temp = alpha*T*zhat+alpha/2*ones(2*par.MT,1);
               
               
               
               
               
               
               
               
               
        end %FIM SWITCH OPCION
               
               
               
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