% author: Xiaojia Lu
% Centre for Wireless Communications (CWC)

function r=MyQuan_LR(symb,N,T,M_QAM)

alpha=2;

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