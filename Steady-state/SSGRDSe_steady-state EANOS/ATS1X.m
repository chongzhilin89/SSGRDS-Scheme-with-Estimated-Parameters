function [EANOS]=ATS1X(m,n,n1,n2,L,L1,L2,L3,dmin,dmax)

% Compute the steady-state ANOS1 of the Case-U SSGRDS Xbar Chart
% using Method 2 as in Khoo et al. (2015) steady-state results
% L3 is from CRL sub-chart

r=sqrt(n1+n2);
c=r/sqrt(n2);
N=4*L3+1;

h=15;
[u,wu]=quadhermite(h);
a=m*(n-1)/2;
b=2/(m*(n-1));
sv=sqrt((a-exp(2*(gammaln(a+0.5)-gammaln(a))))*b);
vmin=max(0,1-5*sv);
vmax=1+5*sv;
[v,wv]=quadlegendre(h,vmin,vmax);
U=ones(h,1)*u';
V=v*ones(1,h);
W=(ones(h,1)*wu').*(wv*ones(1,h));
FU=normpdf(U);
FV=2*V.*gampdf(V.^2,m*(n-1)/2,2/(m*(n-1)));
sn1=sqrt(n1);
sn2=sqrt(n2);
sn1mn=sqrt(n1/(m*n));
sn2mn=sqrt(n2/(m*n));
sn12=sqrt(n1+n2);


[xi,we]=lgwt_table(9); % Abscissas and Weight of Legendre-Gauss Quadrature
EANOSsum=0;
for ii=1:9;
delta=((dmax-dmin)/2)*xi(ii)+((dmax+dmin)/2);

Pr1=normcdf(U*sn1mn+V*L1-delta*sn1)-normcdf(U*sn1mn-V*L1-delta*sn1);
Pr1_0=normcdf(U*sn1mn+V*L1)-normcdf(U*sn1mn-V*L1);

[za,wza]=quadlegendre(h,-L,-L1);
[zb,wzb]=quadlegendre(h,L1,L);
for iu=1:h;
  uiu=u(iu);
  for iv=1:h;
    viv=v(iv);
    fza=viv*normpdf(uiu*sn1mn+viv*za-delta*sn1);
    fzb=viv*normpdf(uiu*sn1mn+viv*zb-delta*sn1);
    p21=sum(wza.*(normcdf(uiu*sn2mn+viv*(L2*sn12-za*sn1)/sn2-delta*sn2)-normcdf(uiu*sn2mn-viv*(L2*sn12+za*sn1)/sn2-delta*sn2)).*fza);
    p22=sum(wzb.*(normcdf(uiu*sn2mn+viv*(L2*sn12-zb*sn1)/sn2-delta*sn2)-normcdf(uiu*sn2mn-viv*(L2*sn12+zb*sn1)/sn2-delta*sn2)).*fzb);
    Pal=sum(wzb.*(normcdf(uiu*sn2mn+viv*(L2*sn12-zb*sn1)/sn2-delta*sn2)).*fzb);
    P1=Pr1(iv,iu);
    P2=p21+p22;
    P=1-P1-P2;
    Q=1-P;

    alfa=(1-normcdf(U(iv,iu)*sn1mn+V(iv,iu)*L1-delta*sn1)-Pal)./P;

    R=sparse(N,N);           
    if (L3==1);
    R(:,2)=Q; R(2,3)=P; R(3,4)=alfa*P;R(5,4)=alfa*P;
    R(3:4,5)=(1-alfa)*P;
    else  
    for i=1:N-1;
       for j=2:N;
          if j==i+1; R(i,j)=Q;
       end;
     end;
    end;
               
    R(L3+1,L3+2)=P;R(L3+1,L3+1)=Q;R(2*L3+1,L3+1)=Q;

    R(3*L3+1,L3+1)=Q;R(4*L3+1,L3+1)=Q;
                
    for i=L3+2:2*L3+1;  
       R(i,2*L3+2)=alfa*P;
    end;
                
    for i=3*L3+2:N;  
       R(i,2*L3+2)=alfa*P;
    end;                
                
    for i=L3+2:3*L3+1;  
        R(i,3*L3+2)=(1-alfa)*P;
     end;
    end;

    fza_0=viv*normpdf(uiu*sn1mn+viv*za);
    fzb_0=viv*normpdf(uiu*sn1mn+viv*zb);
    p21_0=sum(wza.*(normcdf(uiu*sn2mn+viv*(L2*sn12-za*sn1)/sn2)-normcdf(uiu*sn2mn-viv*(L2*sn12+za*sn1)/sn2)).*fza_0);
    p22_0=sum(wzb.*(normcdf(uiu*sn2mn+viv*(L2*sn12-zb*sn1)/sn2)-normcdf(uiu*sn2mn-viv*(L2*sn12+zb*sn1)/sn2)).*fzb_0);
    Pal_0=sum(wzb.*(normcdf(uiu*sn2mn+viv*(L2*sn12-zb*sn1)/sn2)).*fzb_0);
    P1_0=Pr1_0(iv,iu);
    P2_0=p21_0+p22_0;
    P_0=1-P1_0-P2_0;
    Q_0=1-P_0;

    alfa_0=(1-normcdf(U(iv,iu)*sn1mn+V(iv,iu)*L1)-Pal_0)./P_0;

    RR=sparse(N,N); 

    if (L3==1);
    RR(:,2)=Q_0; RR(2,3)=P_0; RR(3,4)=alfa_0*P_0;RR(5,4)=alfa_0*P_0;
    RR(3:4,5)=(1-alfa_0)*P_0;
    else  
    for i=1:N-1;
       for j=2:N;
          if j==i+1; RR(i,j)=Q_0;
       end;
     end;
    end;
               
    RR(L3+1,L3+2)=P_0;RR(L3+1,L3+1)=Q_0;RR(2*L3+1,L3+1)=Q_0;

    RR(3*L3+1,L3+1)=Q_0;RR(4*L3+1,L3+1)=Q_0;
                
    for i=L3+2:2*L3+1;  
       RR(i,2*L3+2)=alfa_0*P_0;
    end;
                
    for i=3*L3+2:N;  
       RR(i,2*L3+2)=alfa_0*P_0;
    end;                
                
    for i=L3+2:3*L3+1;  
        RR(i,3*L3+2)=(1-alfa_0)*P_0;
     end;
    end;


    PP=RR./(sum(RR,2)*ones(1,N));         
    PR=PP'-speye(N);
    PR(1,:)=1;        
    uu=sparse(N,1);
    uu(1)=1;
    q=PR\uu;
    ARL_calc= sum((speye(N)-R)'\q);

    P3=normcdf(U(iv,iu)*sn1mn-V(iv,iu)*L1-delta*sn1)-normcdf(U(iv,iu)*sn1mn-V(iv,iu)*L-delta*sn1)+normcdf(U(iv,iu)*sn1mn+V(iv,iu)*L-delta*sn1)-normcdf(U(iv,iu)*sn1mn+V(iv,iu)*L1-delta*sn1);

    ANOS1P(iv,iu)=ARL_calc*(n1+n2*P3);
  end;
end;

ANOS=sum(sum(W.*FU.*FV.*ANOS1P));

EANOS=we(ii)*(1/(dmax-dmin))*ANOS;
 EANOSsum=EANOSsum+EANOS;
end;
EANOS=((dmax-dmin)/2)*EANOSsum;



