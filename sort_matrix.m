function [P,Hordenada,Q,R] = sort_matrix(H,opcion) %P->matriz Permutação, H->matriz
% Input
%    matrix : H matrix to be sorted
%    opcion : Pega as diferentes opções
% Output
%    P : matriz Permutação
%    H : matriz Saída
%    Hordenada = H*P -> Sae ordenada por ordem ascendente em suas normas


switch opcion 
     
    case 1, %Hordenada = H*P -> É ordenada por ordem ascendente em suas normas cada vetor coluna de H
         P = eye(size(H));
         norma = sqrt(sum(H.^2));
         [nor,index] = sortrows(norma);
         Hordenada = H(:,index);
         P = P(:,index);
         [Q,R] = qr(Hordenada);
         
    case 2, %SQRD Algoritmo-> |R(Nt,Nt)|>....>|R(1,1)|
         R = zeros(size(H)); Q = H; P = eye(size(H));
         norma = sqrt(sum(Q.^2));
         normaux = norma;
         for ii=1:size(H,2)
             [valor,ki] = min(normaux);
             R(:,[ii ki]) = R(:,[ki ii]);
             P(:,[ii ki]) = P(:,[ki ii]);
             Q(:,[ii ki]) = Q(:,[ki ii]);
             R(ii,ii) = norm(Q(:,ii));
             Q(:,ii) = Q(:,ii)/R(ii,ii);
             for j=ii+1:size(H,2)
                 R(ii,j) = Q(:,ii)'*Q(:,j);
                 Q(:,j) = Q(:,j) - R(ii,j)*Q(:,ii);
             end
             normaux = sqrt(sum(Q.^2));     
             normaux(1:ii) = Inf;
         end
         Hordenada = Q*R(1:size(H,2),:);    %Notas que P não foi incluida
         R = R(1:size(H,2),:);
         P = P(1:size(H,2),:);
end %fim Switch