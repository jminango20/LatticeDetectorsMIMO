% author: Shahriar Shahabuddin
% Centre for Wireless Communications (CWC)

function r=MyQuan_CLR(symb_in,N,T_CLR,M_QAM)

symb = [real(symb_in);imag(symb_in)];
T = [real(T_CLR),-imag(T_CLR);imag(T_CLR),real(T_CLR)]; 

alpha=2;
% [M,N]=size(h_lr_temp);
%V = [ -(sqrt(M_QAM)-1) : 2 : (sqrt(M_QAM)-1) ];

z=round(1/alpha*symb-1/2*pinv(T)*ones(N,1));

r_temp=alpha*T*z+alpha/2*ones(N,1);

switch M_QAM
    case 2
        V=[1;-1];
         

    case 4
        V=[1;-1];
        
    case 8
        V=[1;-1;3;-3];
       
        
    case 16
        V=[1;-1;3;-3];
        
        
    case 32
        V=[1;-1;3;-3;5;-5];       
       
        
    case 64
        V=[1;-1;3;-3;5;-5;7;-7]; 
        
end
        L=length(V);
        for kk=1:length(r_temp)
            for ii=1:L
                temp(ii)=(r_temp(kk)-V(ii))^2;                
            end
            indx=find(min(temp)==temp);
            r(kk,:)=V(indx(1));
        end