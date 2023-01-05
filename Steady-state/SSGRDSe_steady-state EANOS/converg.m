function ok=converg(ATS0, tau)

% check whether the ATS0 meets the requirement

if abs(ATS0-tau)<0.01
    ok=1;
else
    ok=0;
end

