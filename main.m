clear; clc; close all;

%% Settings
SeparateLyapunov = 1;               % 1 (0) for using separate (common) Lyapunov functions for
SwLogic = 1;                        % 0: non-switching, 1: arbitrary switching, 2:ADT switching 
BiSearch = 1;                      
MaxNum_bs = 10;                     
muSearchRange = [10^0.001 10^2];    

% for ADT switching
tau = 1000;
mu = 2;
lambda = 1-mu^(-1/tau);
% tau = - log(mu)/log(1-lambda);


if SwLogic ~=2
    BiSearch = 0;
end

% test different partitioning
% Theta1 =[-0.7 1];
% Theta1 =[-0.7 0.2;0.1  1];
Theta1 =[-0.7 0;-0.1 0.5;0.4 1];
% Theta1 =[-0.7 0;-0.1  1];

% Theta1 = [-1 1];
% Theta1 = [-1 0.1;0 1];
% Theta1 = [-1 -0.2;-0.3 0.5;0.4 1];

solvername       = 'sedumi';    
%% Parameter set
Theta1 = [-0.7 0; -0.1 0.5; 0.4 1];
Theta2 = 0;
delta  = 1;
thetaMin = Theta1(1,1);
thetaMax = Theta1(end,end);

if SwLogic == 0
    Theta1 = [Theta1(1,1) Theta1(end,end)];
end

%% Plant
syms theta1 theta2
Gasym.A   = 0.5*[1-theta1  0  -2+theta1;
                 2-theta1 -1   1-theta1;
                -1+theta1  1-3*theta1 -theta1];
Gasym.B1  = [0; 1-theta1; theta1];
Gasym.B2  = [1; 0; 0];
Gasym.C1  = [1 1 1];
Gasym.D11 = 0;
Gasym.D12 = 0;
Gasym.C2  = eye(3);
Gasym.D21 = 0;
Gasym.D22 = 0;

F_theta   = @(x) [1 x(1)];
FthetaNum = [1 1];

%% Region info
regnum1 = size(Theta1,1);
regnum2 = max(size(Theta2,1),1);
regnum  = regnum1*regnum2;

REGID = zeros(regnum,3);
for regid = 1:regnum
    regid1 = mod(regid-1,regnum1)+1;
    regid2 = ceil(regid/regnum1);
    REGID(regid,:) = [regid regid1 regid2];
end

%% Admissible regions
for r1 = 1:regnum1
    [ThetaT{r1}, ThetaPlusT{r1}] = AdmRegDiscrete(Theta1(r1,:), thetaMin, thetaMax, delta);
end
if sum(Theta2(:)) ~= 0
    for r2 = 1:regnum2
        [Theta2T{r2}, Theta2PlusT{r2}] = AdmRegDiscrete(Theta2(r2,:), thetaMin, thetaMax, delta);
    end
else
    Theta2T = {0}; Theta2PlusT = {0};
end

%% Solve
[Gamopt, gamVec, Kopt] = SLPV_SF_discrete_yalmip(Gasym, Theta1, Theta2, ThetaT, ThetaPlusT, Theta2T, Theta2PlusT, REGID, F_theta, FthetaNum, SwLogic, mu, lambda, SeparateLyapunov, solvername);


if ~exist('Kopt','var') || isempty(Kopt)
    error('Kopt does not exist or is empty. Run main.m first.');
end

%% -------------------- Basic parameters --------------------
regnum   = size(Kopt.G,3);
n        = size(Kopt.G,1);          % state dimension
nu       = size(Kopt.L,1);          % number of control inputs
nF       = size(Kopt.L,3);          % number of basis functions
thetaMin = Theta1(1,1);
thetaMax = Theta1(end,end);



Ngrid      = 300;
theta_plot = linspace(thetaMin, thetaMax, Ngrid);
K_all      = zeros(nu, n, Ngrid, regnum);   % K(:,:,k,r)
X_eigs     = zeros(n, Ngrid, regnum);
Xhat_eigs  = zeros(n, Ngrid, regnum);

for r = 1:regnum
    Ginv = inv(Kopt.G(:,:,r));
    for k = 1:Ngrid
        th  = theta_plot(k);
        Ft  = F_theta(th);          % [1, theta]
        
        % L(theta)
        Lth = zeros(nu,n);
        for j = 1:nF
            Lth = Lth + Ft(j)*Kopt.L(:,:,j,r);
        end
        K_all(:,:,k,r) = Lth * Ginv;
        
        % X(theta) and Xhat(theta)
        Xth  = zeros(n);
        Xhth = zeros(n);
        for j = 1:nF
            Xth  = Xth  + Ft(j)*Kopt.X(:,:,j,r);
            Xhth = Xhth + Ft(j)*Kopt.Xhat(:,:,j,r);
        end
        X_eigs(:,k,r)    = sort(real(eig(Xth)),'descend');
        Xhat_eigs(:,k,r) = sort(real(eig(Xhth)),'descend');
    end
end


%% -------------------- Time-domain simulation --------------------
Nsim = 80;                         
t    = 0:Nsim-1;
% scenarios for theta(k)
theta_traj = zeros(3, Nsim);
theta_traj(1,:) = thetaMin + (thetaMax-thetaMin)*(0.5+0.5*sin(2*pi*t/40)); % scenarios 1: slow sinusoid
theta_traj(2,:) = thetaMin + (thetaMax-thetaMin)*(t/Nsim);                  % scenarios 2: ramp
theta_traj(3,:) = thetaMin + (thetaMax-thetaMin)*rand(1,Nsim);              % scenarios 3: random

x0 = [1; -0.5; 0.8];          
figure('Name','Time-domain Simulation','Color','w','Position',[30 30 1250 850]);
for scen = 1:3
    x     = zeros(n, Nsim);
    u     = zeros(nu,Nsim);
    reg   = zeros(1,Nsim);
    x(:,1)= x0;
    
    for k = 1:Nsim-1
        th = theta_traj(scen,k);
        % determine active region (simple, no hysteresis for visualisation)
        r = find(th >= Theta1(:,1) & th <= Theta1(:,2), 1, 'first');
        if isempty(r), r = regnum; end
        reg(k) = r;
        
        Ft  = F_theta(th);
        Lth = zeros(nu,n);
        for j = 1:nF
            Lth = Lth + Ft(j)*Kopt.L(:,:,j,r);
        end
        Kth = Lth / Kopt.G(:,:,r);
        u(:,k) = Kth * x(:,k);
        
        Ga = AugPltEval(Gasym, th);
        x(:,k+1) = Ga.A*x(:,k) + Ga.B2*u(:,k);   % no disturbance for pure free response
    end
    reg(Nsim) = reg(Nsim-1);
    
    % plots for this scenario
    subplot(3,4,(scen-1)*4+1);
    plot(t, theta_traj(scen,:), 'k','LineWidth',1.5); hold on;
    for r = 1:regnum
        idx = reg==r;
        plot(t(idx), theta_traj(scen,idx), '.', 'Color',colors(r,:), 'MarkerSize',8);
    end
    ylabel('\theta(k)'); title(sprintf('Scenario %d – Scheduling',scen));
    grid on; ylim([thetaMin-0.1 thetaMax+0.1]);
    
    subplot(3,4,(scen-1)*4+2);
    plot(t, x', 'LineWidth',1.4); grid on;
    ylabel('x'); title('States');
    legend('x1','x2','x3','Location','best');
    
    subplot(3,4,(scen-1)*4+3);
    plot(t, u', 'LineWidth',1.4); grid on;
    ylabel('u'); title('Control input');
    
    subplot(3,4,(scen-1)*4+4);
    stairs(t, reg, 'LineWidth',1.6); grid on;
    ylabel('Region'); ylim([0.5 regnum+0.5]);
    title('Active region');
    xlabel('time step k');
end
