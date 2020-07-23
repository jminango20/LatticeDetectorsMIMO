%Paper From LR-AidedDetection Towards MLD in MIMO System para 16QAM

function xr = LRdec_improved_16QAM(yr,Hr,par)

xr_hat = zeros(size(Hr,2));

  for k = 1:size(Hr,2)
      Haux = Hr; 
      Haux(:,k)=[];
      E = zeros(size(Hr,1),size(Hr,2)-1); %É uma matriz de dimensão 2K x (2K-1)      
      
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
      % simbolo -3            
      y1 = yr - Hr(:,k)*(-3);
      xr_hat(:,2*k-1) = E * LRdec( y1 , Haux , par );
      xr_hat(k,2*k-1) = -3;
      % simbolo -1      
      y2 = yr - Hr(:,k)*(-1);
      xr_hat(:,2*k) = E * LRdec( y2 , Haux , par );
      xr_hat(k,2*k) = -1;
      % simbolo +1      
      y3 = yr - Hr(:,k)*(1);
      xr_hat(:,4*k-1) = E * LRdec( y3 , Haux , par );
      xr_hat(k,4*k-1) = 1;
      % simbolo +3
      y4 = yr - Hr(:,k)*(3);
      xr_hat(:,4*k) = E * LRdec( y4 , Haux , par );
      xr_hat(k,4*k) = 3;

      
  end
  
  Hrxr = Hr*xr_hat;
  yr = repmat(yr,1,size(Hrxr,2));
  xr_distancia = sum(abs(yr-Hrxr).^2);
  [valor,pos]=min(xr_distancia);
  xr = xr_hat(:,pos);