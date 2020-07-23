%   Modified LLL (MLLL) Lattice Reduction Algorithm 
%   Original Paper: A Customized Lattice Reduction Multiprocessor for MIMO Detection
%   by Shahriar Shahabuddin et. al.
%
%   T = MLLL(H,Tx);
%
%   Input:
%   - H          Complex-valued channel matrix
%   - Tx         Number of Transmit Antenna
%   Output:
%   - T          Unimodular matrix needed for LR aided detectors
%
%   by Shahriar Shahabuddin
%   Email: shahriar.cwc@gmail.com

function T = MLLL(H,Tx)

[Q,R] = qr(H); 

%initialize
lambda = 3/4;
mm = Tx;
T     = eye(mm);                        
S_flag = 0;
iter = 0;
% Set the number of iterations
iter_max = 5;

while S_flag == 0 && iter<iter_max

    S_flag = 1;
    iter = iter+1;

  %size reduction
for k = 2:mm 
   for l=k-1:-1:1
      mu = round(R(l,k)/R(l,l));      
      if abs(mu) ~= 0
         R(1:l,k) = R(1:l,k) - mu * R(1:l,l);
         T(:,k)   = T(:,k)   - mu * T(:,l);
      end
   end
   
    temp_div = norm(R(k-1:k,k));
   %Siegel Condition
   if lambda*abs(R(k-1,k-1)) > abs(R(k,k))
      
      %exchange columns k and k-1     
        temp_r=R(:,k-1);
        temp_t=T(:,k-1);
        
        R(:,k-1)=R(:,k);
        R(:,k)=temp_r;
        T(:,k-1)=T(:,k);
        T(:,k)=temp_t;

      % Given's rotation
      alpha     = R(k-1,k-1) / temp_div;        
      beta     = R(k,k-1)   / temp_div;
      Theta = [alpha' beta; -beta alpha];
      
      R(k-1:k,k-1:end) = Theta * R(k-1:k,k-1:end);
      S_flag = 0;
   end
      k = k+1;
end
end