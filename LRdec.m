function xr = LRdec(yr,Hr,par) %Paper From LR-AidedDetection Towards MLD in MIMO System
 Nt = size(Hr,2); 
 P = LLL(Hr,Nt/2);
 Hred = Hr*P;
 xr_hat = pinv(Hred)*yr;
 %xr_hat = P* round(( pinv(Hred)*yr - pinv(P) * (ones(Nt,1)) )) + (ones(Nt,1)); %Está no Paper
 xr = MyQuan_LR(xr_hat,Nt,P,length(par.symbols)); %Chamo a função MyQuan_LR
 
 
 
