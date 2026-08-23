function Gplt = AugPltEval(Gplt_sym, Par)
nPar = numel(Par);
if nPar >= 1
    theta1 = Par(1); 
end
if nPar >= 2
    theta2 = Par(2); 
end
if nPar >= 3
    theta3 = Par(3); 
end

Gplt.A   = double(subs(Gplt_sym.A));
Gplt.B1  = double(subs(Gplt_sym.B1));
Gplt.B2  = double(subs(Gplt_sym.B2));
Gplt.C1  = double(subs(Gplt_sym.C1));
Gplt.D11 = double(subs(Gplt_sym.D11));
Gplt.D12 = double(subs(Gplt_sym.D12));
Gplt.C2  = double(subs(Gplt_sym.C2));
Gplt.D21 = double(subs(Gplt_sym.D21));
Gplt.D22 = double(subs(Gplt_sym.D22));
end