%   Lenstra-Lenstra-Lovasz (LLL) Lattice Reduction Algorithm
%   Original Paper: Near-maximum-likelihood detection of MIMO systems
%   using MMSE-Based Lattice-Reduction, by Dirk Wubben et al.
%
%   T = LLL(H_real,Tx);
%
%   Input:
%   - H_real     Real valued channel matrix
%   - Tx         Number of Transmit Antenna
%   Output:
%   - T          Unimodular matrix needed for LR aided detectors
%
%   by Shahriar Shahabuddin
%   Email: shahriar.cwc@gmail.com


function T = LLL(H_real,Tx)

[Q,R]=qr(H_real);

% Real domain operation, thats why multiplied by 2.
mm=2*Tx;

%initialize
R_lr=R; T=eye(mm);
% lamda = 3/4 was proposed in the original LLL paper.
lamda=3/4;

k=2;

% size reduction operation
while k<=mm,
    for l=k-1:-1:1
        mu=round( R_lr(l,k)/R_lr(l,l) );
        if mu ~=0
            R_lr(1:l,k)=R_lr(1:l,k)-mu*R_lr(1:l,l);
            T(:,k)=T(:,k)-mu*T(:,l);
        end
    end

% LLL reduction

    if lamda*R_lr(k-1,k-1)^2>R_lr(k,k)^2+R_lr(k-1,k)^2%&&k<=mm
        temp_r=R_lr(:,k-1);
        temp_t=T(:,k-1);%exchange columns k and k-1
        
        R_lr(:,k-1)=R_lr(:,k);
        R_lr(:,k)=temp_r;
        T(:,k-1)=T(:,k);
        T(:,k)=temp_t;
        
        %Givens rotation R_lr(k,k-1) = 0
        alpha=R_lr(k-1,k-1)/norm(R_lr(k-1:k,k-1));
        beta=R_lr(k,k-1)/norm(R_lr(k-1:k,k-1));
        Theta=[alpha,beta;-beta,alpha];
    
        R_lr(k-1:k,k-1:mm)=Theta*R_lr(k-1:k,k-1:mm);
        k=max(k-1,2);
    else
        k=k+1;
    end
    
end