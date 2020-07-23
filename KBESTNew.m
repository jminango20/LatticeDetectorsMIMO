function xhat_quan = KBESTNew(H_real,yreal,par,s_real,w,sigma)  %POSSO TIRAR NA ULTIMA PARTE QUE APROXIMO AO BOUNDARY A FIM DE VER COMO MELHORE o DESEMPENHO LEVANDO AO BOUNDARY
K = 2;  %K -> Ramos Melhores ; N -> Vizinhos ; onde N>K
%Adaptação Na Lattice
alpha = 2; Lattice = ([real(par.symbols);imag(par.symbols)] + ones(2,par.M))/alpha;
%Simbolos Transmitidos no dominio da Lattice
x = ( s_real + ones(2*par.MT,1) )/alpha; 
H_real = [ H_real ; sigma*eye(2*par.MT) ];  %H_real Extendida
T = LLL(H_real,par.MT); 
z = pinv(T)*x;
H_tiu = H_real*T;

%% POSSO ESCOLHER QUALQUER DOS DOIS FORMAS DE DESCOMPOSIÇÂO DE QR
[Q,R] = qr(H_tiu); P = eye(size(H_tiu)); %As duas Linhas vão juntas quando eu usar elas--ATENÇÂO.
R = R(1:size(H_tiu,2),:);
P = P(1:size(H_tiu,2),:);
Q = Q(1:2*(par.MT+par.MR),1:2*par.MT);
dR = 2*(diag(R) >= 0) - 1; % avoid using sign which can return 0, destroying qr-Diag(R) é positiva
n = length(dR);
Q(:,1:n) = bsxfun(@times, Q(:,1:n), dR'); R(1:n,:) = bsxfun(@times, R(1:n,:), dR);

%ORDENADO
%[P,H_tiu,Q,R] = sort_matrix(H_tiu,2);  %Chamo ao sort_matrix que é SQRD
%Para a qr do Matlab


%% Vetor Simbolos Transmitido
zordenado = P*z;
y = H_tiu*zordenado + [(w)/alpha ; -sigma*zordenado] ;
y = Q'*y;

zk = nan( 2*par.MT+1 , K );
cost = zeros( 2*par.MT+1 , K);
len = 1;

prod = 0; jj=1;
   
for n = 2*par.MT:-1:1
%       %[z(1:len , n ) , cost(1:len , n ) ] = FindK( z(1:len , n+1 ) , cost(1:len , n+1 ) );
%       
    for i=1:len
          
           for j=n+1:2*par.MT
               prod(jj) = R(n,j)*zk(j,i); %LINHA 2 , atenção
               jj = jj+1;
           end
           eliminar = sum(prod); prod=0;         
           ri(i) = y(n) - eliminar;
           zi(i) = round(ri(i)/R(n,n));
           child(:,i) = [zi(i) ; zk(n+1,i) ];
           childcost(i) = cost(n+1,i) + (ri(i) - R(n,n)*zi(i))^2;
           step(i) = sign(ri(i)/R(n,n) - zi(i));
    end %fim for i       
    q = childcost;  %Ver q=childcost(i) onde i=1 até len - q as the keys
      
    for k=1:K
        [a,idx] = min(q);   %Encontro o indice i associado à minimo key in q
        zk(n,k) = child(1,idx); %NOTA
        cost(n,k) = childcost(idx);
        zi(i) = zi(idx) + step(i);
        child(:,i) = [zi(i) ; zk(n+1,i) ];
%       child = [ child(:,i) [zi(i);z(n+1,i)] ];
        childcost(i) = cost(n+1,i) + (ri(i) - R(n,n)*zi(i))^2;
        step(i) = -step(i) - sign(step(i));
        q = [q childcost(i)];
    end
    len = K;
end %fim For n
zhat = zk(1:2*par.MT,:);

yy = repmat(y,1,size(zhat,2));
[minimo,posicion] = min(sum(abs(yy - R*zhat).^2));
zhat_minimo = zhat(:,posicion);    
shat = alpha*T*zhat_minimo - ones(2*par.MT,1);
r_temp = shat;   %Levar nos limites do Alfabeto usado
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

