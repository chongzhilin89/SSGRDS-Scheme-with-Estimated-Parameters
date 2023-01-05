
% Optimization program for Case-U Side Sensitive Group Run DS chart(Steady state)
% L3 is from CRL sub-chart

tic;
m=1000000   %change as required
n=3        %change as required
dmin = 0.2 %change as required
dmax = 1.4 %change as required
tau=370     %change as required

ATS1minopt=inf;
for n1 = 1:n-1;  %for1
 for n2 = 1:2*n; %for2
 if (n1+n2)>n; %if1
    ATS1min=inf;
    ARL_0=tau/n;
    alfa_0=1/ARL_0;
    L1min=icdf('Normal',1-(n-n1)/(2*n2)-alfa_0/2,0,1);
    L1 = L1min;
    L=icdf('Normal',(n-n1)/(2*n2)+normcdf(L1,0,1),0,1);

    L3ok=1;
    if n < 6
    L3=4; 
    else 
        L3=13;
    end
    while L3ok==1 %control the L3 loop

    [L2,ATS0]=tsearch_L2(m,n,n1,n2,L,L1,L3,tau);
    ATS1=ATS1X(m,n,n1,n2,L,L1,L2,L3,dmin,dmax);

    if abs(ATS1-ATS1min)<0.005; % Terminate the loop if no significant reduction in ATS1
    L3ok=0;
    end;
            
     if ATS1 < ATS1min  %if2(control the while loop)                   
        ATS1min=ATS1;
        ATS0min=ATS0;
        n1ot=n1; n2ot=n2;
        Lot=L; L1ot=L1; L2ot=L2; L3ot=L3;
        L3=L3+1;
        Optimal_solution=[L3ot,n1ot,n2ot,Lot,L1ot,L2ot,ATS0min,ATS1min]
  
     else
        L3ok=0;
     end; %end if2          
    end; %end while    

    if ATS1min<ATS1minopt %if3                 
       n1ot1=n1ot; n2ot1=n2ot;
       Lot1=Lot; L1ot1=L1ot; L2ot1=L2ot; L3ot1=L3ot;
       ATS0minopt=ATS0min; ATS1minopt=ATS1min;
    end; %end if3           
  end; %end if1
 end; %end for2
end; %end for1

m
n
dmin
dmax
L3ot1
FinalOpt_n1_n2_L_L1_L2=[n1ot1,n2ot1,Lot1,L1ot1,L2ot1] % final optimal solution 
ATS0minopt
ATS1minopt

toc;
