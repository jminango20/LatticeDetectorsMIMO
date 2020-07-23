%r = -5.87; a = r - floor(r); b = 1 - a; N = 20; c = 3.16;

function xhat = PDF(r,c)

     a = r - floor(r);
     b = 1 - a;
     N = 10;  %Default
     for i=0:N-1
         s(i+1) = exp(-c*((a+i)^2)) + exp(-c*((b+i)^2));    
     end
     s_linha = sum(s);
     
     posicion = 1;
     %Para exp com a
     for ii=1:(N)   
         P(posicion) = exp( -c * ( ( a + N - ii )^2 ) )/s_linha;
         q(posicion) = floor(r) - N + ii;
         posicion = posicion + 1;        
     end     
     %Para exp com b
     for ii=N:-1:1    
         P(posicion) = exp( -c * ( ( b + N - ii )^2 ) )/s_linha;
         q(posicion) = floor(r) + N - ii + 1;
         posicion = posicion + 1;
     end 
     P;
     [v,pos]=max(P);
     xhat = q(pos);



%OUTRO
%a = r - floor(r);
%b = 1-a;
%for i=0:200
%    sparcial(i+1) = exp(-c*(a+i)^2) + exp(-c*(b+i)^2);
%end
%s = sum(sparcial);
%for i=0:200  %Pus 200
%    qa(i+1) = floor(r) - i;
%    Pa(i+1) = (exp(-c*(a+i)^2))/s;
%    qb(i+1) = floor(r) + 1 + i;
%    Pb(i+1) = (exp(-c*(b+i)^2))/s;
%end
%P = [Pa Pb];
%q = [qa qb];
%[v,pos] = max(P);
%xhat = q(pos);

