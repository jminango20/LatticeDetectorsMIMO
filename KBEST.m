function [xhat_quan,posicion] = KBEST(H_real,yreal,par,s_real,w,sigma,opcion)  %POSSO TIRAR NA ULTIMA PARTE QUE APROXIMO AO BOUNDARY A FIM DE VER COMO MELHORE o DESEMPENHO LEVANDO AO BOUNDARY

K = 4; N = 4;  %K -> Ramos Melhores ; N -> Vizinhos ; onde N>K
%Adaptação Na Lattice
alpha = 2; 
Lattice = ([real(par.symbols);imag(par.symbols)] + ones(2,par.M))/alpha;

%Simbolos Transmitidos no dominio da Lattice
x = ( s_real + ones(2*par.MT,1) )/alpha; 
if opcion == 1
    T = LLL(H_real,par.MT); 
else
    H_real = [ H_real ; sigma*eye(2*par.MT) ];  %H_real Extendida
    T = LLL(H_real,par.MT); 
end
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
if opcion == 1
    y =  H_tiu*zordenado + w/alpha;  %Na Lattice
else
    y = H_tiu*zordenado + [(w)/alpha ; -sigma*zordenado] ;
end
y = Q'*y;
%Onde vou salvar os K best
C = NaN(2*par.MT,K); 
Dbest = zeros(2*par.MT,K);

%% Algoritmo
%Aplicação do SIC na Layer 2*par.MT
l = 2*par.MT;
zhat(l,1) = y(l)/R(l,l);
zhatquan(l,1) = round(zhat(l,1));
%Vizinhos = zhatquan(l,1); deucl = abs(Vizinhos - zhat(l,1))^2; %Quantização
Vizinhos = zhatquan(l,1); 
deucl = abs(y(l) - R(l,l)*Vizinhos)^2;
for ii=1:N
    Vizinhos(2*ii) = zhatquan(l,1) - ii;    
    %deucl(2*ii) = abs(Vizinhos(2*ii) - zhat(l,1))^2;
    deucl(2*ii) = abs(y(l) - R(l,l)*Vizinhos(2*ii))^2;
    Vizinhos(2*ii+1) = zhatquan(l,1) + ii;    
    %deucl(2*ii+1) = abs(Vizinhos(2*ii+1) - zhat(l,1))^2;
    deucl(2*ii+1) = abs(y(l) - R(l,l)*Vizinhos(2*ii+1))^2;
end
[valor,idx] = sort(deucl);
C(l,:) = Vizinhos(idx(1:K));  %K melhores vetores
Dbest(l,:) = deucl(idx(1:K));  %distância euclideana dos K melhores vetores
d = Dbest;

%Aqui começa o codigo para as (2*par.MT Layer - 1)
for k=l-1:-1:1 %k-> são os Layer 
    for i=1:K  %i->Candidato parcial
        zhat(k,i) = (y(k) - R(k,k+1:2*par.MT)*C(k+1:2*par.MT,i)) / R(k,k) ;       
        zhatquan(k,i) = round(zhat(k,i));
        Vizinhos(i,1) = zhatquan(k,i);
        NVetores(:,1,i) = [Vizinhos(i,1);C(k+1:2*par.MT,i)];        
        %dparcial(i,1) = d(k+1,i) + abs(Vizinhos(i,1) - zhat(k,i)).^2;
        dparcial(i,1) = d(k+1,i) + abs(y(k) - R(k,k)*Vizinhos(i,1) - R(k,k+1:2*par.MT)*C(k+1:2*par.MT,i) ).^2;
        for ii=1:N  %Gero os N pontos vizinhos
            Vizinhos(i,2*ii) = zhatquan(k,i) - ii;
            NVetores(:,2*ii,i) = [Vizinhos(i,2*ii);C(k+1:2*par.MT,i)];
            dparcial(i,2*ii) = d(k+1,i) + abs(y(k) - R(k,k)*Vizinhos(i,2*ii) - R(k,k+1:2*par.MT)*C(k+1:2*par.MT,i)).^2;            
            Vizinhos(i,2*ii+1) = zhatquan(k,i) + ii;    
            NVetores(:,2*ii+1,i) = [Vizinhos(i,2*ii+1);C(k+1:2*par.MT,i)];
            dparcial(i,2*ii+1) = d(k+1,i) + abs(y(k) - R(k,k)*Vizinhos(i,2*ii) - R(k,k+1:2*par.MT)*C(k+1:2*par.MT,i)).^2;                        
        end
    end %fim K best
    %Pegar os Kbest melhores com menor norma Euclideana
    for kk = 1:K
        [r,c] = find(dparcial==min(min(dparcial)));
        C(k:2*par.MT,kk) = NVetores(:,c(1),r(1));  %Atenção
        Dbest(k,kk) = dparcial(r(1),c(1));
        d(k,kk) = dparcial(r(1),c(1));
        dparcial(r(1),c(1)) = inf;  
    end
    Vizinhos=0;NVetores=[];dparcial=0;   
end %fim de Layer

%% Detecção
zhat = C;   %Todos os K-best simbolos vectores
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
