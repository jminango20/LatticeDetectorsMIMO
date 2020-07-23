%Paper From LR-AidedDetection Towards MLD in MIMO System para 16QAM - com
%Coset Decomposition

function xr = LRdec_improved_16QAM_coset(yr,Hr,par)

xr_hat = zeros(size(Hr,2));

  for k = 1:size(Hr,2)
      Haux = Hr; 
      %Haux(:,k)=[];
      %E = zeros(size(Hr,1),size(Hr,2)-1); %É uma matriz de dimensão 2K x (2K-1)            
      %for ii=1:size(E,1) 
      %    for jj=1:size(E,2)
      %        if ii<k && jj==ii
      %            E(ii,jj)=1;
      %        elseif ii>k && jj==(ii-1)
      %            E(ii,jj)=1;
      %        else
      %            E(ii,jj)=0;
      %        end
      %    end
      %end
      
      v = ones(size(Hr,2),1); v(k)= 2; %O que faz o 2-> no programa
      % xr E -1            
      b = -1;
      y1 = yr - Hr(:,k)*(-1);
      xr_hat(:,2*k-1) = LRdec( y1 , Haux*diag(v) , par );
      xr_hat(k,2*k-1) = xr_hat(k,2*k-1) + b;
      % xr E +1
      b = 1;      
      y2 = yr - Hr(:,k)*(1);
      xr_hat(:,2*k) = LRdec( y2 , Haux*diag(v) , par );
      xr_hat(k,2*k) = xr_hat(k,2*k) + b;

      
  end
  
  Hrxr = Hr*xr_hat;
  yr = repmat(yr,1,size(Hrxr,2));
  xr_distancia = sum(abs(yr-Hrxr).^2);
  [valor,pos]=min(xr_distancia);
  xr = xr_hat(:,pos);