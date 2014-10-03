function alpha=galoistable(i,n)
for i=1:3
    alpha{1,i}=dec2bin(i-1,n);
end
alpha{1,4}=dec2bin(4,n);
alpha{1,5}=dec2bin(8,n);
alpha{1,6}=dec2bin(3,n);
alpha{1,7}=dec2bin(6,n);
alpha{1,8}=dec2bin(12,n);
alpha{1,9}=dec2bin(11,n);
alpha{1,10}=dec2bin(5,n);
alpha{1,11}=dec2bin(10,n);
alpha{1,12}=dec2bin(7,n);
alpha{1,13}=dec2bin(14,n);
alpha{1,14}=dec2bin(15,n);
alpha{1,15}=dec2bin(13,n);
alpha{1,16}=dec2bin(9,n);

zero=dec2bin(0,n)
one=dec2bin(1,n)

k=3
for j=1:14
    alphapower(1,j)=alpha(1,k)
    k=k+1;
end



%%% To calculate PA=2*G,GIVEN Xp=1,Yp=alpha.^6 %%%%%%%
XP=1;G=[1,alphapower(1,6)]
LAMBDA1=dec2bin(bitxor(bin2dec(one),bin2dec(alphapower(1,6))),4)
for j=1:14
    if(LAMBDA1==alphapower{1,j})
       LAMBDA1 = j
        break
    else
        j=j+1;
    end
end
newj=j*2
if(newj>=15)
    newj=mod(newj,15)
end
XR=dec2bin(bitxor(bin2dec(alphapower(1,13)),bin2dec(alphapower(1,newj))),4);
XR=dec2bin(bitxor(bin2dec(XR),bin2dec(alphapower(1,4))),4);
Xr1=bin2dec(XR)

YR=dec2bin(bitxor(bin2dec(alphapower(1,13)),bin2dec(one)),4);
YR=bin2dec(YR);
YR=Xr1*0;
Yr1=XP.^2+YR
PA{1,1}=Xr1;
PA{1,2}=Yr1;
fprintf('PA co ordinates are:alphapower%d and alphapower%d \n',Xr1,Yr1)
%%% To calculate PB=3*G,GIVEN Xp=1,Yp=alpha.^6 %%%%%%%
NR=dec2bin(bitxor(bin2dec(alphapower(1,6)),bin2dec(one)),4)
NR1=bin2dec(NR)
DR1=dec2bin(bitxor(bin2dec(one),bin2dec(zero)),4)
DR1=bin2dec(DR1)
if(DR1==1)
LAMBDA2=NR1
else if(DR1==0)
LAMBDA2=1/0
    else
    LAMBDA2=NR1+(15-DR1)
    end
end
LAMBDA2=dec2bin(LAMBDA2)
for m=1:14
    if(LAMBDA2==alphapower{1,m})
       LAMBDA2 = m
        break
    else
        m=m+1;
    end
end
XR2=dec2bin(bitxor(bin2dec(alphapower(1,13)),bin2dec(alphapower(1,newj))),4);
XR2=dec2bin(bitxor(bin2dec(XR2),bin2dec(one)),4);
XR2=dec2bin(bitxor(bin2dec(XR2),bin2dec(zero)),4);
XR2=dec2bin(bitxor(bin2dec(XR2),bin2dec(alphapower(1,4))),4);
Xr2=bin2dec(XR2)
YR2=LAMBDA2
PB{1,1}=Xr2
PB{1,2}=YR2
fprintf('THE co-ordinate  points are:alphapower%d and alphapower%d \n',Xr2,YR2)
%%% To calculate nAPB=2*Pb,GIVEN Xp=1,Yp=alphapower13 %%%%%%%
NR2=dec2bin(bitxor(bin2dec(alphapower(1,YR2)),bin2dec(XR2)),4)
NR2=bin2dec(NR2)

DR2=bin2dec(XR2)
if(DR2==1)
LAMBDA3=NR2
else if(DR2==0)
LAMBDA3=1/0;
    else
    LAMBDA3=NR2*(15-DR2);
    end
end
LAMBDA3=dec2bin(LAMBDA3)
for l=1:14
    if(LAMBDA3==alphapower{1,l})
       LAMBDA3 = l
        break
    else
        l=l+1;
    end
end
LAMBDAsq=LAMBDA3.*2

XR3=dec2bin(bitxor(bin2dec(alphapower(1,LAMBDAsq)),bin2dec(alphapower(1,LAMBDA3))),4);
for O=1:14
    if(XR3==alphapower{1,O})
       XR3 = O
        break
    else
        O=O+1;
    end
end
XR3=dec2bin(bitxor(bin2dec(alphapower(1,XR3)),bin2dec(alphapower(1,4))),4)
if (XR3=='0000')
    XR3=0
else if( XR3=='0001')
        XR3=1
    else
for P=1:14
    if(XR3==alphapower{1,P})
       XR3 = P
        break
    else
        P=P+1;
    end
end
    end
end

YR3=dec2bin(bitxor(bin2dec(alphapower(1,LAMBDA3)),bin2dec(one)),4);
if (YR3=='0000')
    YR3=0
else if( YR3=='0001')
        YR3=1
    else
for Q=1:14
    if(YR3==alphapower{1,Q})
       YR3 = Q
        break
    else
        Q=Q+1;
    end
end
    end
end
YR3=YR3*XR3;
NEWXp=Xr2.^2;
YR3=dec2bin(bitxor(YR3,NEWXp),4)
if (YR3=='0000')
    YR3=0;
else if( YR3=='0001')
        YR3=1;
    else
for R=1:14
    if(YR3==alphapower{1,R})
       YR3 = R;
        break
    else
        R=R+1;
    end
end
    end
end
Ka=[XR3 YR3]
fprintf('SECRET KEY Ka points are:alphapower%d and alphapower%d \n',XR3,YR3)
% To calculate Pm+na.Pb %
Pm{1,1}=alphapower{1,10}
Pm{1,2}=alphapower{1,8}
NUM4=dec2bin(bitxor(bin2dec(Pm{1,2}),YR3),4)
if(NUM4==0)
    NUM4=0;
else if(NUM4==1)
        NUM4=1;
    else
for S=1:14
    if(NUM4==alphapower{1,S})
       NUM4 = S
        break
    else
        S=S+1;
    end
end
    end
end
  DEN4=dec2bin(bitxor(bin2dec(Pm{1,1}),XR3),4)
  if(DEN4==0)
    DEN4=0;
else if(DEN4==1)
        DEN4=1;
    else
for S=1:14
    if(DEN4==alphapower{1,S})
       DEN4 = S
        break
    else
        S=S+1;
    end
end
    end
  end
  
if(DEN4==1)
LAMBDA4=NUM4
else if(DEN4==0)
LAMBDA4=1/0;
    else
    LAMBDA4=NUM4+(15-DEN4)
    end
end
if LAMBDA4>=15
    LAMBDA4=mod(LAMBDA4,15)
    LAMBDA4SQ=LAMBDA4*2
else
LAMBDA4SQ=LAMBDA4*2
end
if(LAMBDA4SQ>=15)
    LAMBDA4SQ=mod(LAMBDA4SQ,15)
end
XR4=dec2bin(bitxor(bin2dec(alphapower(1,LAMBDA4SQ)),bin2dec(alphapower(1,LAMBDA4))),4);
if(XR4==0)
    XR4=0
else if(XR4==1)
        XR4=1
    else
for S=1:14
    if(XR4==alphapower{1,S})
       XR4 = S;
        break
    else
        S=S+1;
    end
end
    end
end
XR4=dec2bin(bitxor(bin2dec(alphapower(1,XR4)),bin2dec(Pm{1,1})),4)
if(XR4==0)
    XR4=0
else if(XR4==1)
        XR4=1
    else
for S=1:14
    if(XR4==alphapower{1,S})
       XR4 = S
        break
    else
        S=S+1;
    end
end
    end
end

 if(XR3==0 || XR3==1)
 XR4a=dec2bin(bitxor(bin2dec(alphapower(1,4)),XR3),4);
 
 else
    XR4a=dec2bin(bitxor(bin2dec(alphapower(1,4)),bin2dec(alphapower(1,XR3))),4);
 end
 if(XR4a==0)
    XR4a=0
else if(XR4a==1)
        XR4a=1
    else
 for T=1:14
    if(XR4a==alphapower{1,T})
       XR4a = T;
        break
    else
        T=T+1;
    end
end
    end
 end
 XR4=dec2bin(bitxor(bin2dec(alphapower(1,XR4a)),bin2dec(alphapower(1,XR4))),4)
 if(XR4==0)
    XR4=0
else if(XR4==1)
        XR4=1
    else
for S=1:14
    if(XR4==alphapower{1,S})
       XR4 = S
        break
    else
        S=S+1;
    end
end
    end
 end

YR4a=dec2bin(bitxor(bin2dec(Pm{1,1}),bin2dec(alphapower(1,XR4))),4);
if(YR4a==0)
    YR4a=0
else if(YR4a==1)
        YR4a=1
    else
for S=1:14
    if(YR4a==alphapower{1,S})
       YR4a = S
        break
    else
        S=S+1;
    end
end
    end
end
YR4b=dec2bin(bitxor(bin2dec(Pm{1,2}),bin2dec(alphapower(1,XR4))),4);
if(YR4b==0)
    YR4b=0
else if(YR4b==1)
        YR4b=1
    else
for S=1:14
    if(YR4b==alphapower{1,S})
       YR4b = S
        break
    else
        S=S+1;
    end
end
    end
end
 YR4c=YR4a+LAMBDA4
 if(YR4c>=15)
   YR4c=mod(YR4c,15)
 end
YR4=dec2bin(bitxor(bin2dec(alphapower(1,YR4c)),bin2dec(alphapower(1,YR4b))),4)
if(YR4==0)
    YR4=0
else if(YR4==1)
        YR4=1
    else
for S=1:14
    if(YR4==alphapower{1,S})
       YR4 = S
        break
    else
        S=S+1;
    end
end
    end
end
%Display the cipher text %
 Cm{1,1}=XR3
 Cm{1,2}=YR3
 Cm{1,3}=XR4
 Cm{1,4}=YR4
  Cm
  fprintf('THE CIPHERED points are:alphapower%d alphapower%d,alphapower%d alphapower%d \n',XR3,YR3,XR4,YR4)

 %-----------------------END OF ENCRYPTION------------------%
 
 
 % DECRYPTION %
 
 if (Xr1==0 && Yr1==1) 
     y2=[Xr1 Yr1]
 end
 % generalise it later %%%%
%    y2=dec2bin(bitxor(bin2dec(alphapower(1,Yr1)),bin2dec(alphapower(1,Xr1))),4);
% for S=1:14
%     if(y2==alphapower{1,S})
%        y2 = S;
%         break
%     else
%         S=S+1;
%     end
% end
%  end
%  y2=dec2bin(bitxor(bin2dec(alphapower(1,y2)),bin2dec(alphapower(1,))),4);
% for S=1:14
%     if(YR4==alphapower{1,S})
%        YR4 = S;
%         break
%     else
%         S=S+1;
%     end
% end

Yr1=Yr1+Xr1
y2=[Xr1 Yr1]
%
if Yr1==1 || Yr1==0
NUM5=dec2bin(bitxor(Yr1,bin2dec(alphapower(1,YR4))),4)
else
 NUM5=dec2bin(bitxor(bin2dec(alphapower(1,Yr1)),bin2dec(alphapower(1,YR4))),4) 
end
for S=1:14
    if(NUM5==alphapower{1,S})
       NUM5 = S
        break
    else
        S=S+1;
    end
end
if Xr1==1 || Xr1==0
DEN5=dec2bin(bitxor(Xr1,bin2dec(alphapower(1,XR4))),4)
else
 DEN5=dec2bin(bitxor(bin2dec(alphapower(1,Xr1)),bin2dec(alphapower(1,XR4))),4) 
end

for S=1:14
    if(DEN5==alphapower{1,S})
       DEN5 = S
        break
    else
        S=S+1;
    end
end
if(DEN5==1)
LAMBDA5=NUM5
else if(DEN5==0)
LAMBDA5=disp('INFINITY');
    else
    LAMBDA5=NUM5+(15-DEN5)
    end
end
LAMBDA5SQ=LAMBDA5*2
if(LAMBDA5SQ>=15)
    LAMBDA5SQ=mod(LAMBDA5SQ,15)
end
XR5=dec2bin(bitxor(bin2dec(alphapower(1,LAMBDA5SQ)),bin2dec(alphapower(1,LAMBDA5))),4);
for S=1:14
    if(XR5==alphapower{1,S})
       XR5 = S;
        break
    else
        S=S+1;
    end
end

XR5=dec2bin(bitxor(bin2dec(alphapower(1,XR5)),bin2dec(alphapower(1,XR4))),4)
for S=1:14
    if(XR5==alphapower{1,S})
       XR5 = S
        break
    else
        S=S+1;
    end
end
 if(Xr1==0 || Xr1==1)
 XR5a=dec2bin(bitxor(bin2dec(alphapower(1,4)),Xr1),4);
 
 else
    XR5a=dec2bin(bitxor(bin2dec(alphapower(1,4)),bin2dec(alphapower(1,Xr1))),4);
 end
 for T=1:14
    if(XR5a==alphapower{1,T})
       XR5a = T;
        break
    else
        T=T+1;
    end
end

 XR5=dec2bin(bitxor(bin2dec(alphapower(1,XR5a)),bin2dec(alphapower(1,XR5))),4)
for S=1:14
    if(XR5==alphapower{1,S})
       XR5 = S
        break
    else
        S=S+1;
    end
end

YR5a=dec2bin(bitxor(bin2dec(alphapower(1,XR5)),bin2dec(alphapower(1,XR4))),4);
for S=1:14
    if(YR5a==alphapower{1,S})
       YR5a = S
        break
    else
        S=S+1;
    end
end
YR5b=dec2bin(bitxor(bin2dec(alphapower(1,XR5)),bin2dec(alphapower(1,YR4))),4)
if(YR5b =='0000')
    YR5b=0
else
   if(YR5b =='0001')
    YR5b=1
    
else
for S=1:14
    if(YR5b==alphapower{1,S})
       YR5b = S
        break
    else
        S=S+1;
    end
end
   end
end
YR5c=YR5a+LAMBDA5;
 if(YR5c>=15)
   YR5c=mod(YR5c,15)
 end
 if(YR5b==1 || YR5b==0)
YR5=dec2bin(bitxor(bin2dec(alphapower(1,YR5c)),YR5b),4)
 
for S=1:14
    if(YR5==alphapower{1,S})
       YR5 = S
        break
    else
        S=S+1;
    end
end
 else 
     
     YR5=dec2bin(bitxor(bin2dec(alphapower(1,YR5c)),bin2dec(alphapower(1,YR5b))),4)
for S=1:14
    if(YR5==alphapower{1,S})
       YR5 = S
        break
    else
        S=S+1;
    end
end
 end
PMD=[XR5 YR5]

fprintf('THE DECIPHERED points are:alphapower%d and alphapower%d \n',XR5,YR5)

