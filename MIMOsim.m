function MIMOsim(varargin)
  % -- Parametros por default/Parametros do Usuario
  if isempty(varargin)
    disp('Usando Parametros Simula��o Por Default...')        
    % Parametros por default
    par.simName = 'ERR_4x4_16QAM'; % Nome da Simula��o (usado para salvar resultados)
    par.runId = 0; % simula��o ID (usado para reproduzir resultados)
    par.MR = 4; % Antenas Receptoras
    par.MT = 4; % Antenas Transmissoras (n�o maiores que MR!) 
    par.mod = 'QPSK'; % Tipo Modula��o: 'BPSK','QPSK','16QAM','64QAM'
    par.SNRdB_list = 0:3:21; % Valores da SNR [dB] a ser simulados
    %par.detector = {'LLL_ZF','LLL_ZF_SIC','LLL_ZF_Fix','ZF','LLL_MMSE','ELR_ZF','ELR_ZF2','Amostragem','LatticeAumento','ImprovedLatticeAumento','ImprovedLatticeAumentoMMSE','KMelhores'}; % Defini��o detectore(s) a serem simulados
    %par.detector = {'LLL_ZF','ZF','LLL_MMSE','ELR_ZF','Amostragem','LatticeAumento','ImprovedLatticeAumento','ImprovedLatticeAumentoMMSE'}; % Defini��o detectore(s) a serem simulados
    %par.detector = {'ELR_ZF2','MODProcedure'}; %N�O SERVEM REVISAR-- Defini��o detectore(s) a serem simulados
    par.detector = {'LRdecimproved2','LRdecimproved','ImprovedLatticeAumentoMMSE','LLL_ZF','ImprovedLatticeAumento'}; 
    
  else      
    disp('usando parametros j� estabelecidos...')    
    par = varargin{1}; % s� argumento de estrutura
  end

  % -- Inizializa��o
  % use runId random seed (enables reproducibility)
  %rng(par.runId); 

  % Alfabeto da Constela��o com Mapeamento Gray-mapped 
  switch (par.mod)
    case 'BPSK',
      par.symbols = [ -1 1 ];
    case 'QPSK', 
      par.symbols = [ -1-1i,-1+1i, ...
                      +1-1i,+1+1i ];
    case '16QAM',
      par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
                      -1-3i,-1-1i,-1+3i,-1+1i, ...
                      +3-3i,+3-1i,+3+3i,+3+1i, ...
                      +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
      par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
                      -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
                      -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
                      -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
                      +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
                      +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
                      +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
                      +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];                         
  end

  % Obten��o average de Energia de Simbolo
  par.Es = mean(abs(par.symbols).^2);   
  % precomputo bits
  par.Q = log2(length(par.symbols)); % numero de bits por simbolo
  par.M = length(par.symbols); % number of simbolos da modula��o
  par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');
  par.n_iterations = 10;     % numero de diferentes canais de H
  par.n_symbols = 500;    % numero de vetores Transmitidos por cada canal de H


  %% Parametros da Simula��o
  n = sqrt(0.5)*(randn(par.MR,par.n_symbols,par.n_iterations)+1i*randn(par.MR,par.n_symbols,par.n_iterations)); %Vetor Ruidos
  H = sqrt(0.5)*(randn(par.MR,par.MT,par.n_iterations)+1i*randn(par.MR,par.MT,par.n_iterations)); %Vetor Canal Rayleigh
  bits = randi([0 1],par.MT,par.Q,par.n_symbols,par.n_iterations);  %Bits Gerados para a Transmiss�o

  %% Parametros Dos Resultados
  res.SER = zeros(length(par.SNRdB_list),length(par.detector)); % symbol error rate
  res.BER = zeros(length(par.SNRdB_list),length(par.detector)); % bit error rate
  res.time = zeros(length(par.detector),1); % Tempo Ejecu��o de Cada Algoritmo

  %% -- COME�A A SIMULA��O  
  %parfor i=1:n_iterations      
  for i=1:par.n_iterations
      i
      
      for d=1:length(par.detector)         
        switch (par.detector{d}) % Escolha Algoritmos de Detec��o

%% DETECTOR MLD            
            case 'MLD', % MLD Detec��o - Tradicional            
                tic;
                [BER , SER] = MLD(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;
                
%% DETECTOR ESFERICO                            
            case 'SD_Studer', % MLD Detec��o - Esferico Studer: Mais Rapido que do Acima            
                tic;
                [BER , SER] = SD_Studer(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;

%% DETECTOR ZERO FORCING            
            case 'ZF', % ZF Detec��o - Tradicional            
                tic;
                [BER , SER] = ZF(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                

            case 'ZF_SIC', % ZF-SIC Detec��o - Tradicional            
                tic;
                [BER , SER] = ZF_SIC(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                

%% DETECTOR MMMSE                            
            case 'MMSE', % MMSE Detec��o - Tradicional            
                tic;
                [BER , SER] = MMSE(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                
                
            case 'MMSE_SIC', % MMSE-SIC Detec��o - Tradicional            
                tic;
                [BER , SER] = MMSE_SIC(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                
                
%% DETECTOR ZF com Redu��o Lattice LLL_ZF                            
            case 'LLL_ZF', % Lattice com ZF Detec��o - Tradicional            
                tic;
                [BER , SER] = LLL_ZF(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    
                
            case 'LLL_ZF_SIC', % Lattice com ZF com SIC Detec��o - Tradicional            
                tic;
                [BER , SER] = LLL_ZF_SIC(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    
                
            case 'LLL_ZF_Fix', % Lattice com ZF com Fix Detec��o - Tradicional            
                tic;
                [BER , SER] = LLL_ZF_Fix(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    

%% DETECTOR MMSE com Redu��o Lattice LLL_MMSE                            
            case 'LLL_MMSE', % Lattice com MMSE Detec��o - Tradicional            
                tic;
                [BER , SER] = LLL_MMSE(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% DETECTOR Elementary Lattice Reduction com ZF
            case 'ELR_ZF', % E_Lattice_Reduction com ZF Detec��o
                tic;
                [BER , SER] = ELR_ZF(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                    

%% DETECTOR Elementary Lattice Reduction 2 com ZF -- REVISAS NO CONVERGE
            case 'ELR_ZF2', % % E_Lattice_Reduction2 com ZF Detec��o 
                tic;
                [BER , SER] = ELR_ZF2(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                    
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% NOVOS DETECTORES                
            case 'MODProcedure',  %Novo Decoder           
                tic;
                [BER , SER] = MODProcedure(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;                    
                
%% NOVOS DETECTORES                
            case 'Amostragem',  %Novo Decoder           
                tic;
                [BER , SER] = Amostragem(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    

%% NOVOS DETECTORES LatticeAumento                
            case 'LatticeAumento',  %Novo Decoder           
                tic;
                [BER , SER] = LatticeAumento(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    

%% NOVOS DETECTORES ImprovedLatticeAumento               
            case 'ImprovedLatticeAumento',  %Novo Decoder           
                tic;
                [BER , SER] = ImprovedLatticeAumento(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    
                
%% NOVOS DETECTORES ImprovedLatticeAumento               
            case 'ImprovedLatticeAumentoMMSE',  %Novo Decoder           
                tic;
                [BER , SER] = ImprovedLatticeAumentoMMSE1(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    
                
%% NOVOS DETECTORES KMELHORES
            case 'KMelhores',  %Novo Decoder           
                tic;
                [BER , SER] = KMelhores(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    
                
%% NOVOS DETECTORES LRImproved - Com Um Elemento a Mudan�a
            case 'LRdecimproved',  %Novo Decoder           
                tic;
                [BER , SER] = LRdecimproved(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    

%% NOVOS DETECTORES LRImproved - Com Dois Elementos a Mudan�a
            case 'LRdecimproved2',  %Novo Decoder           
                tic;
                [BER , SER] = LRdecimproved2(par.n_symbols , 1 , bits(:,:,:,i) , n(:,:,i) , H(:,:,i) , par);
                res.SER(:,d) = res.SER(:,d) + SER/par.n_iterations; 
                res.BER(:,d) = res.BER(:,d) + BER/par.n_iterations; 
                res.time(d,1) = res.time(d,1) + toc;    
                
                
            otherwise,            
                error('tipo de par.detector n�o definido.')      
        end
      end %fim d->detector
      
  end %fim i->n_iteration
  
  %% MOSTRA DOS RESULTADOS
  save([ par.simName '_' num2str(par.runId) ],'par','res');    
  
%% -- Mostra Resultados (Genera��o Matlab plot)  -> BER
  marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
  figure(1)
  for d=1:length(par.detector)
      if d==1
          semilogy(par.SNRdB_list,res.BER(:,d),marker_style{d},'LineWidth',2)
          hold on      
      else          
          semilogy(par.SNRdB_list,res.BER(:,d),marker_style{d},'LineWidth',2)
      end    
  end
  %hold off
  grid on
  xlabel('SNR [dB]','FontSize',12)
  ylabel('bit error rate (BER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-5 1])
  legend(par.detector,'FontSize',12)
  title(['Compara��o de Varios Detectores em ',num2str(par.MT),'x',num2str(par.MR),' sistema, ', par.mod, ' simbolos']);
  set(gca,'FontSize',12)

  %% -- Mostra Resultados (Genera��o Matlab plot)  -> SER
  figure(2)
  for d=1:length(par.detector)
      if d==1
          semilogy(par.SNRdB_list,res.SER(:,d),marker_style{d},'LineWidth',2)
          hold on      
      else          
          semilogy(par.SNRdB_list,res.SER(:,d),marker_style{d},'LineWidth',2)
      end    
  end
  %hold off
  grid on
  xlabel('SNR [dB]','FontSize',12)
  ylabel('Symbol Error Rate (SER)','FontSize',12)
  axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-5 1])
  legend(par.detector,'FontSize',12)
  title(['Compara��o de Varios Detectores em ',num2str(par.MT),'x',num2str(par.MR),' sistema, ', par.mod, ' simbolos']);
  set(gca,'FontSize',12)
  
  %% -- Mostra Resultados (Genera��o Matlab plot)  -> Tempo Simula��o
  figure(3)
  bar(res.time/60)
  set(gca,'XTickLabel',par.detector);
  xlabel('Detectores'), ylabel('Tempo Computacional (min)')
