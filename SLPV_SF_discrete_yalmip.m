function [Gamopt, gamVec, Kopt] = SLPV_SF_discrete_yalmip(Gasym, Theta1, Theta2, ...
    ThetaT, ThetaPlusT, Theta2T, Theta2PlusT, REGID, F_theta, FthetaNum, ...
    SwLogic, mu, lambda, MultiLyapunov, solvername)

% Discrete-time Switching LPV State-Feedback with Separate Lyapunov
% YALMIP version (equivalent to the original LMI Lab code)

if nargin < 15
    solvername = 'mosek';   % default
end

if isempty(Theta2T) || (iscell(Theta2T) && sum(Theta2T{1})==0)
    GSParaNum = 1;
else
    GSParaNum = 2;
end

% Plant dimensions
if GSParaNum == 1
    theta0 = ThetaT{1}(1);
else
    theta0 = [ThetaT{1}(1); Theta2T{1}(1)];
end
Ga = AugPltEval(Gasym, theta0);
n  = size(Ga.A,1);
nw = size(Ga.B1,2);
nu = size(Ga.B2,2);
nz = size(Ga.C1,1);

regnum  = max(REGID(:,1));
regnum1 = max(REGID(:,2));
regnum2 = max(REGID(:,3));

% ========== Decision variables ==========
Gam = sdpvar(1);

for r = 1:regnum
    gam(r) = sdpvar(1);
    G{r}   = sdpvar(n, n, 'full');
    
    for j = 1:sum(FthetaNum)
        X{j,r}    = sdpvar(n, n);
        if MultiLyapunov
            Xhat{j,r} = sdpvar(n, n);
        else
            Xhat{j,r} = X{j,r};
        end
        L{j,r} = sdpvar(nu, n, 'full');
    end
end

Constraints = [];
Objective   = Gam;

% ========== LMIs for each region ==========
for r = 1:regnum
    rid1 = REGID(r,2);
    rid2 = REGID(r,3);
    
    thT   = ThetaT{rid1};
    thpT  = ThetaPlusT{rid1};
    th2T  = Theta2T{rid2};
    th2pT = Theta2PlusT{rid2};
    
    for i1 = 1:length(thT)
        theta1  = thT(i1);
        theta1p = thpT(i1);
        
        for i2 = 1:length(th2T)
            theta2  = th2T(i2);
            theta2p = th2pT(i2);
            
            theta  = [theta1; theta2];
            thetap = [theta1p; theta2p];
            
            Ga = AugPltEval(Gasym, theta);
            A  = Ga.A;  B1 = Ga.B1;  B2 = Ga.B2;
            C1 = Ga.C1; D11= Ga.D11; D12= Ga.D12;
            
            Ft  = F_theta(theta);
            Ftp = F_theta(thetap);
            
            % Build X(theta), X(theta+), L(theta)
            Xth  = zeros(n);
            Xthp = zeros(n);
            Lth  = zeros(nu,n);
            Xhth = zeros(n);
            Xhthp= zeros(n);
            
            for j = 1:sum(FthetaNum)
                Xth   = Xth   + Ft(j)  * X{j,r};
                Xthp  = Xthp  + Ftp(j) * X{j,r};
                Lth   = Lth   + Ft(j)  * L{j,r};
                Xhth  = Xhth  + Ft(j)  * Xhat{j,r};
                Xhthp = Xhthp + Ftp(j) * Xhat{j,r};
            end
            
            % ----- Local performance LMI (dilated discrete) -----
            % Equivalent to H_analysis / H_synthesis in the paper
            vu = 1;   % for performance
            if SwLogic == 2 && ~MultiLyapunov
                vu = 1 - lambda;
            end
            
            M11 = vu*(G{r}+G{r}') - vu*Xth;
            M21 = A*G{r} + B2*Lth;
            M22 = Xthp;
            M31 = C1*G{r} + D12*Lth;
            M33 = gam(r)*eye(nz);
            M42 = B1';
            M43 = D11';
            M44 = gam(r)*eye(nw);
            
            LMI_perf = [M11 , M21' , M31' , zeros(n,nw);
                M21 , M22  , zeros(n,nz) , M42';
                M31 , zeros(nz,n) , M33 , M43';
                zeros(nw,n) , M42 , M43 , M44];
            
            Constraints = [Constraints, LMI_perf >= 0];
            
            % ----- Stability LMI with Xhat (when separate) -----
            if MultiLyapunov && SwLogic ~= 0
                vu_g = 1;
                if SwLogic == 2
                    vu_g = 1 - lambda;
                end
                
                M11s = vu_g*(G{r}+G{r}') - vu_g*Xhth;
                M21s = A*G{r} + B2*Lth;
                M22s = Xhthp;
                
                LMI_stab = [M11s , M21s';
                    M21s , M22s];
                Constraints = [Constraints, LMI_stab >= 0];
            end
        end
    end
    
    % gam(r) <= Gam
    Constraints = [Constraints, gam(r) <= Gam, gam(r) >= 1e-6];
end

% ========== Switching surface conditions ==========
if SwLogic ~= 0
    for r2 = 1:regnum2
        RegId1 = sort(REGID(REGID(:,3)==r2, 2));
        RegId  = (r2-1)*regnum1 + RegId1;
        
        for ii = 1:length(RegId1)-1
            for th2 = unique(Theta2(r2,:))
                % direction 1
                th1 = Theta1(RegId1(ii), end);
                Ft  = F_theta([th1; th2]);
                Xa = zeros(n); Xb = zeros(n);
                for j = 1:sum(FthetaNum)
                    Xa = Xa + Ft(j)*Xhat{j, RegId(ii)};
                    Xb = Xb + Ft(j)*Xhat{j, RegId(ii+1)};
                end
                Constraints = [Constraints, Xa <= mu * Xb];
                
                % direction 2
                th1 = Theta1(RegId1(ii+1), 1);
                Ft  = F_theta([th1; th2]);
                Xa = zeros(n); Xb = zeros(n);
                for j = 1:sum(FthetaNum)
                    Xa = Xa + Ft(j)*Xhat{j, RegId(ii+1)};
                    Xb = Xb + Ft(j)*Xhat{j, RegId(ii)};
                end
                Constraints = [Constraints, Xa <= mu * Xb];
            end
        end
    end
end

% ========== Solve ==========
ops = sdpsettings('solver', solvername, 'verbose', 1, ...
    'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS', 1e-8, ...
    'sedumi.eps', 1e-8);

sol = optimize(Constraints, Objective, ops);

if sol.problem == 0
    Gamopt = value(Gam);
    gamVec = zeros(1,regnum);
    for r = 1:regnum
        gamVec(r) = value(gam(r));
        Kopt.G(:,:,r) = value(G{r});
        for j = 1:sum(FthetaNum)
            Kopt.L(:,:,j,r)    = value(L{j,r});
            Kopt.X(:,:,j,r)    = value(X{j,r});
            Kopt.Xhat(:,:,j,r) = value(Xhat{j,r});
        end
    end
else
    Gamopt = inf;
    gamVec = inf(1,regnum);
    Kopt   = [];
    warning('Solver returned problem code: %d', sol.problem);
end
end