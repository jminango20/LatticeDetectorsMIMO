%Paper From LR-AidedDetection Towards MLD in MIMO System para 4QAM

function xr = LRdec_improved(par,s_real,H_real,w)

opcion = 2;
valor = 1;   %Salva os parametros de Cond.Number e Ortogonalidade

    if opcion ~= 2
        yr = H_real*s_real + w;
    else
        yr = H_real*s_real + (w*(0.9));
    end
    xr_hat = zeros(size(H_real,2));


    for k = 1:size(H_real,2)
        Haux = H_real; 
        Haux(:,k)=[];
        E = zeros(size(H_real,1),size(H_real,2)-1); %É uma matriz de dimensão 2K x (2K-1)      
        for ii=1:size(E,1) 
            for jj=1:size(E,2)
                if ii<k && jj==ii                
                    E(ii,jj)=1;                
                elseif ii>k && jj==(ii-1)
                    E(ii,jj)=1;
                else                    
                    E(ii,jj)=0;
                end
            end
        end        
        yred(:,1) = yr - H_real(:,k)*(1);
        yred(:,2) = yr + H_real(:,k)*(1);
        %% LLL
        T = LLL(Haux , size(Haux,2)/2);
        Hred = Haux*T;
        Cond_Haux(valor) = cond(Haux);
        Cond_Hred(valor) = cond(Hred);
        OrtDeficiencia_Haux(valor) = 1 - ( det(Haux'*Haux) / prod((sum(Haux.^2))) );
        OrtDeficiencia_Hred(valor) = 1 - ( det(Hred'*Hred) / prod((sum(Hred.^2))) );
        valor = valor + 1;
        
        %% Criação Possiveis Vetores Solução
        xr_hat(:,2*k-1) = E * MyQuan_LR(pinv(Hred)*yred(:,1),size(Haux,2),T,length(par.symbols)); %Chamo a função MyQuan_LR
        xr_hat(k,2*k-1) = 1;

        xr_hat(:,2*k) = E * MyQuan_LR(pinv(Hred)*yred(:,2),size(Haux,2),T,length(par.symbols)); %Chamo a função MyQuan_LR
        xr_hat(k,2*k) = -1;
    end
        
    switch opcion
        
        case 1  %Distância Euclideana 
            Hrxr = H_real*xr_hat;
            yr = repmat(yr,1,size(Hrxr,2));
            xr_distancia = sum(abs(yr-Hrxr).^2);
            [dist_menor,pos]=min(xr_distancia);
            xr = xr_hat(:,pos);       
        
        case 2  %Distância Euclideana e Estresada
            Hrxr = H_real*xr_hat;
            yr = repmat(yr,1,size(Hrxr,2));
            xr_distancia = sum(abs(yr-Hrxr).^2);
            [dist_menor,pos]=min(xr_distancia);
            xr = xr_hat(:,pos);       
            %Estressada a xr. -> Vou mudando elemento por elemento de xr e avalio
            %sua distância.
            for ii=1:size(xr,1)        
                if xr(ii,1)==1
                   xr(ii,1) = -1;
                   dist_aux = sum(abs( yr(:,1) - H_real*xr ).^2);
                   if dist_aux > dist_menor
                      xr(ii,1) = 1;
                   else
                      dist_menor = dist_aux;
                   end
                elseif xr(ii,1)==-1                
                    xr(ii,1) = 1;
                    dist_aux = sum(abs( yr(:,1) - H_real*xr ).^2);
                    if dist_aux > dist_menor
                        xr(ii,1) = -1;
                    else                        
                        dist_menor = dist_aux;
                    end                    
                end                
            end %fim FOR
        
        case 3  %Análisis de Signo em xrhat e Estressada
            for jj=1:size(xr_hat,1)
                uno_positivo = sum(xr_hat(jj,:)==1);
                uno_negativo = sum(xr_hat(jj,:)==-1);
                
                if uno_positivo>uno_negativo
                    xr_celso(jj,1) = 1;
                elseif uno_positivo<uno_negativo
                    xr_celso(jj,1) = -1;
                elseif uno_positivo==uno_negativo
                    xr_celso(jj,1) = 1;
                end
            end %fim FOR jj
            dist_menor = sum(abs( yr(:,1) - H_real*xr_celso ).^2);
            
            %Estressada a xr_celso. -> Vou mudando elemento por elemento de xr e avalio
            %sua distância.
            %for ii=1:size(xr_celso,1)        
            %    if xr_celso(ii,1)==1
            %       xr_celso(ii,1) = -1;
            %       dist_aux = sum(abs( yr(:,1) - H_real*xr_celso ).^2);
            %       if dist_aux > dist_menor
            %          xr_celso(ii,1) = 1;
            %       else
            %          dist_menor = dist_aux;
            %       end
            %    elseif xr_celso(ii,1)==-1                
            %        xr_celso(ii,1) = 1;
            %        dist_aux = sum(abs( yr(:,1) - H_real*xr_celso ).^2);
            %        if dist_aux > dist_menor
            %            xr_celso(ii,1) = -1;
            %        else                        
            %            dist_menor = dist_aux;
            %        end                    
            %    end                
            %end %fim FOR
            xr = xr_celso;
            
        
        case 4  %Análisis de Signo em xrhat e Estresada com Toma de Decisão por Probabilidades        
            duda = 0;
            posicion = 1;
            for jj=1:size(xr_hat,1)
                uno_positivo = sum(xr_hat(jj,:)==1);
                Prob_uno_positivo = uno_positivo/size(xr_hat,2);
                uno_negativo = sum(xr_hat(jj,:)==-1);
                Prob_uno_negativo = uno_negativo/size(xr_hat,2);               
                if Prob_uno_positivo >= 0.9
                    xr_duda(jj,1) = 1;
                elseif Prob_uno_negativo >= 0.9
                    xr_duda(jj,1) = -1;
                else 
                    xr_duda(jj,1) = 1;  %Atenção coloco isso só para ter um valor em dist_menor
                    duda(posicion,1) = jj;
                    posicion = posicion+1; 
                end
            end %fim FOR jj
            dist_menor = sum(abs( yr(:,1) - H_real*xr_duda ).^2);
            
            %Estressada a xr_duda. -> Vou mudando elemento por elemento de xr_duda e avalio
            %sua distância.
            if duda ~= 0
                for ii=1:size(duda,1)        
                    if xr_duda(duda(ii),1)==1
                       xr_duda(duda(ii),1) = -1;
                       dist_aux = sum(abs( yr(:,1) - H_real*xr_duda ).^2);
                       if dist_aux > dist_menor
                          xr_duda(duda(ii),1) = 1;
                       else
                          dist_menor = dist_aux;
                       end
                    elseif xr_duda(duda(ii),1)==-1                
                        xr_duda(duda(ii),1) = 1;
                        dist_aux = sum(abs( yr(:,1) - H_real*xr_duda ).^2);
                        if dist_aux > dist_menor
                            xr_duda(duda(ii),1) = -1;
                        else                        
                            dist_menor = dist_aux;
                        end                    
                    end                
                end %fim FOR
            end %fim IF
            xr = xr_duda;
    
    end %FIM SWITCH
        
        
        
        
        
    