%Paper From LR-AidedDetection Towards MLD in MIMO System para 4QAM - AQUI
%PEGO OU ASSUMO DOIS CANDIDATOS

function xr = LRdec_improved2(par,s_real,H_real,w)

opcion = 1;
pos_sol = 1;


    yr = H_real*s_real + w;
    xr_hat = zeros(size(H_real,2));

    
    for k = 1:size(H_real,2)
        for k1 = k+1:size(H_real,2)
            pos = 1;
            Haux = H_real; 
            Haux(:,k)=[]; Haux(:,k1-1)=[];
    
                       
            Eaux = eye(size(H_real,1)-2,size(H_real,2)-2); %É uma matriz de dimensão 2K x (2K-1)          
            E(k,:)=zeros(1,size(H_real,2)-2);
            E(k1,:)=zeros(1,size(H_real,2)-2);
            for jj=1:size(H_real,1)
                if k~=jj
                    if k1~=jj
                        E(jj,:) = Eaux(pos,:); 
                        pos = pos + 1;
                    end
                end
            end            
        
        
            
            yred(:,1) = yr - H_real(:,k)*(1) - H_real(:,k1)*(1) ;
            yred(:,2) = yr - H_real(:,k)*(1) - H_real(:,k1)*(-1) ;
            yred(:,3) = yr - H_real(:,k)*(-1) - H_real(:,k1)*(1) ;
            yred(:,4) = yr - H_real(:,k)*(-1) - H_real(:,k1)*(-1) ;
            
            %LLL
            T = LLL(Haux , size(Haux,2)/2);
            cond(Haux);
            Hred = Haux*T;
            cond(Hred);
            
            %1)
            xr_hat(:,pos_sol) = E*MyQuan_LR(pinv(Hred)*yred(:,1),size(Haux,2),T,length(par.symbols)); %Chamo a função MyQuan_LR
            xr_hat([k k1],pos_sol) = [1;1];
            pos_sol = pos_sol + 1;

            %2)
            xr_hat(:,pos_sol) = E*MyQuan_LR(pinv(Hred)*yred(:,2),size(Haux,2),T,length(par.symbols)); %Chamo a função MyQuan_LR
            xr_hat([k k1],pos_sol) = [1;-1];
            pos_sol = pos_sol + 1;
            
            %3)
            xr_hat(:,pos_sol) = E*MyQuan_LR(pinv(Hred)*yred(:,3),size(Haux,2),T,length(par.symbols)); %Chamo a função MyQuan_LR
            xr_hat([k k1],pos_sol) = [-1;1];
            pos_sol = pos_sol + 1;

            %4)
            xr_hat(:,pos_sol) = E*MyQuan_LR(pinv(Hred)*yred(:,4),size(Haux,2),T,length(par.symbols)); %Chamo a função MyQuan_LR
            xr_hat([k k1],pos_sol) = [-1;-1];
            pos_sol = pos_sol + 1;
            
            
        end %fim k1
    end %fim k
    
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
            for ii=1:size(xr_celso,1)        
                if xr_celso(ii,1)==1
                   xr_celso(ii,1) = -1;
                   dist_aux = sum(abs( yr(:,1) - H_real*xr_celso ).^2);
                   if dist_aux > dist_menor
                      xr_celso(ii,1) = 1;
                   else
                      dist_menor = dist_aux;
                   end
                elseif xr_celso(ii,1)==-1                
                    xr_celso(ii,1) = 1;
                    dist_aux = sum(abs( yr(:,1) - H_real*xr_celso ).^2);
                    if dist_aux > dist_menor
                        xr_celso(ii,1) = -1;
                    else                        
                        dist_menor = dist_aux;
                    end                    
                end                
            end %fim FOR
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
