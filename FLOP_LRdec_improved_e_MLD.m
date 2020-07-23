clc; clear all
NT = 2:10; 
M = 4; %Modulação 4QAM

for ii=1:size(NT,2)
    
    MLDFlop(ii) = ( 8*NT(ii)^2 + 4*NT(ii) - 1) * (M)^(2*NT(ii));
    
    %LRdec_improved(ii) = (2*NT(ii))*( (2*NT(ii)) * (2*NT(ii))^3 * log(2*NT(ii)) ) + ( 8*NT(ii)^2 + 4*NT(ii) - 1)*(2*2*NT(ii));
    LRdec_improved(ii) = (2*NT(ii)) * ( NT(ii) ) * ( (2*NT(ii)) * (2*NT(ii))^3 * log(2*NT(ii)) ) + ( 8*NT(ii)^2 + 4*NT(ii) - 1)*(2*2*NT(ii));
    
end
figure(2)
semilogy(NT,MLDFlop,'bs-');
hold on
semilogy(NT,LRdec_improved,'mv:');
grid on
title('Complexidade')
xlabel('Número Antenas Nr x Nt')
ylabel('Número de Flops')
legend('MLD','Proposto')