function xhat_quan = KBESTNew2(H_real,yreal,par,s_real,w,sigma,opcion) 
    K = 2;  %Numero de possiveis Candidatos

    y = (yreal + H_real*ones(2*par.MT,1))/2;
    T = LLL(H_real,par.MT);
    H = H_real*T;
    [Q,R]=qr(H);
    dR = 2*(diag(R) >= 0) - 1; % avoid using sign which can return 0, destroying qr-Diag(R) é positiva
    n = length(dR);
    Q(:,1:n) = bsxfun(@times, Q(:,1:n), dR'); R(1:n,:) = bsxfun(@times, R(1:n,:), dR);
    z = Q'*y;
    n = 2*par.MT;
    %% III Inicialização 
    yhat(n,1) = z(n)/R(n,n);
    s(n,1) = round(yhat(n,1));
    D(n,1) = (z(n,1) - R(n,n)*s(n,1))^2;
    for ii=1:K/2  %Observar depois o K
        s(n,2*ii) = s(n,1) - ii;    
        D(n,2*ii) = (z(n,1) - R(n,n)*s(n,2*ii))^2;        
        s(n,2*ii+1) = s(n,1) + ii;    
        D(n,2*ii+1) = (z(n,1) - R(n,n)*s(n,2*ii+1))^2;    
    end
    
%% IV Expanssão e Sort    
   for l=n-1:-1:1       
       for ll=1:size(s,2)    
           yhat(l,ll) = (z(l) - R(l,l+1:2*par.MT)*s(l+1:2*par.MT,ll))/R(l,l);  
           C(ll,1) = round(yhat(l,ll));
           Dtemp(ll,1) = D(l+1,ll) + ( z(l) - R(l,l+1:2*par.MT)*s(l+1:2*par.MT,ll) - R(l,l)*C(ll,1) )^2;    
           for ii=1:K  %Observar depois o K
               C(ll,2*ii) = C(ll,1) - ii;    
               Dtemp(ll,2*ii) = D(l+1,ll) + (z(l) - R(l,l+1:2*par.MT)*s(l+1:2*par.MT,ll) - R(l,l)*C(ll,2*ii) )^2;
               C(ll,2*ii+1) = C(ll,1) + ii;    
               Dtemp(ll,2*ii+1) = D(l+1,ll) + (z(l) - R(l,l+1:2*par.MT)*s(l+1:2*par.MT,ll) - R(l,l)*C(ll,2*ii+1) )^2;
           end          
       end %fim ll       
       saux = s;
       for k=1:K+1
           [r,c] = find(Dtemp==min(min(Dtemp)));
           s(l:2*par.MT,k) = [C(r,c);saux(l+1:2*par.MT,r)];
           D(l:2*par.MT,k) = Dtemp(r,c);
           Dtemp(r,c) = inf;
       end %fim K
       Dtemp = 0;       
   end %fim l
   
   zz = repmat(z,1,size(s,2));
   [minimo,posicion] = min(sum(abs(zz - R*s).^2));
   shat_minimo = s(:,posicion);    
   shat = 2*T*shat_minimo - ones(2*par.MT,1);
   %Levar nos limites do Alfabeto usado
   r_temp = shat;   
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
   L = length(V);
   for kk=1:length(r_temp)
       for iii=1:L
           temp(iii)=(r_temp(kk)-V(iii))^2;                
       end       
       indx=find(min(temp)==temp);
       r(kk,:)=V(indx(1));
   end
   xhat_quan = r;
