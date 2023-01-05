function [L2,ATS0]=tsearch_L2(m,n,n1,n2,L,L1,L3,tau)

% search L2 to make ATS0=tau (zero state)

step=0.1;
L2=1.2;
x1=L2;
f1=ATS0X(m,n,n1,n2,L,L1,x1,L3);

ok=converg(f1, tau);
if ok==1
    ATS0=f1;
    L2=x1;
    return
end

if f1 > tau  % decide search direction
    step=-step;
end

% bracket solution
while ok==0
    x2=x1+step;
    f2=ATS0X(m,n,n1,n2,L,L1,x2,L3);
    
    ok=converg(f2, tau);
    if ok==1
        ATS0=f2;
        L2=x2;
        return
    end
    
    if (f1-tau)*(f2-tau) > 0
        step=step.*1.2;
        x1=x2;
        f1=f2;
        ok=0;
    else
        ok=1;
    end
end

% swap x1 with x2 and f1 with f2
if x1 > x2
    tempx=x1;
    x1=x2;
    x2=tempx;
    
    tempf=f1;
    f1=f2;
    f2=tempf;
end

% finalize solution by halving method
while (1)
    %find a candidate value for L by linear interpolating x1 nad x2
    aa=(f2-f1)./(x2-x1);
    bb=f1-(aa.*x1);
    L2=(tau-bb)/aa;
    f=ATS0X(m,n,n1,n2,L,L1,L2,L3);
    
    ok=converg(f, tau);
    if ok==1
        ATS0=f;
        return
    end
    
    % update x1 and x2 by L
    if f > tau
        x2=L2;
        f2=f;
    else
        x1=L2;
        f1=f;
    end
end
    
    