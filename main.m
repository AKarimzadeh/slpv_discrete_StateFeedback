clear; clc; close all;

%% Settings
SeparateLyapunov = 1;               % 1 (0) for using separate (common) Lyapunov functions for
SwLogic = 1;                        % 0: non-switching, 1: arbitrary switching, 2:ADT switching 
BiSearch = 1;                       % do a bisection-search for determining mu value in ADT switching
MaxNum_bs = 10;                     % maximum number of trials for bisection-search 
muSearchRange = [10^0.001 10^2];    % search range for mu in average-dwell-time switching

% for ADT switching
tau = 1000;
mu = 2;
lambda = 1-mu^(-1/tau);
% tau = - log(mu)/log(1-lambda);

% no need for bisearch if not for ADT switching
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

solvername       = 'sedumi';    % 'mosek' | 'sedumi' | 'sdpt3'
%% Parameter set (same as original)
Theta1 = [-0.7 0; -0.1 0.5; 0.4 1];
Theta2 = 0;
delta  = 1;
thetaMin = Theta1(1,1);
thetaMax = Theta1(end,end);

if SwLogic == 0
    Theta1 = [Theta1(1,1) Theta1(end,end)];
end

%% Plant (original discrete example)
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
[Gamopt, gamVec, Kopt] = SLPV_SF_discrete_yalmip(Gasym, Theta1, Theta2, ...
    ThetaT, ThetaPlusT, Theta2T, Theta2PlusT, REGID, F_theta, FthetaNum, ...
    SwLogic, mu, lambda, SeparateLyapunov, solvername);

%% Results
fprintf('\n========== YALMIP Discrete-time SLPV Results ==========\n');
fprintf('Solver: %s\n', solvername);
fprintf('Global gamma = %.6f\n', Gamopt);
fprintf('Local gammas = '); fprintf('%.4f  ', gamVec); fprintf('\n');

%% ========================================================================
%  VISUALIZATION & INTUITIVE ANALYSIS of Switching LPV State-Feedback
%  Run this AFTER main.m has finished successfully
% ========================================================================

if ~exist('Kopt','var') || isempty(Kopt)
    error('Kopt does not exist or is empty. Run main.m first.');
end

%% -------------------- 0. Basic parameters --------------------
regnum   = size(Kopt.G,3);
n        = size(Kopt.G,1);          % state dimension
nu       = size(Kopt.L,1);          % number of control inputs
nF       = size(Kopt.L,3);          % number of basis functions
thetaMin = Theta1(1,1);
thetaMax = Theta1(end,end);

fprintf('\n================================================================\n');
fprintf('           SWITCHING LPV – COMPREHENSIVE RESULTS\n');
fprintf('================================================================\n');
fprintf('Global gamma (Gamopt)          = %.6f\n', Gamopt);
fprintf('Local gammas (gamVec)          = '); fprintf('%.4f  ', gamVec); fprintf('\n');
fprintf('Number of regions              = %d\n', regnum);
fprintf('State dimension n              = %d\n', n);
fprintf('Input dimension nu             = %d\n', nu);
fprintf('Basis functions for L/X        = %d\n', nF);
fprintf('Separate Lyapunov              = %d\n', SeparateLyapunov);
fprintf('Switching logic (0/1/2)        = %d\n', SwLogic);
fprintf('================================================================\n\n');

%% -------------------- 1. Reconstruct K(theta) on a dense grid --------------------
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

%% -------------------- 2. Plot Scheduled Gains K_ij(theta) --------------------
figure('Name','Scheduled Gains K(theta)','Color','w','Position',[50 50 1200 750]);
colors = lines(regnum);
for i = 1:nu
    for j = 1:n
        subplot(nu,n,(i-1)*n+j);
        hold on; grid on; box on;
        for r = 1:regnum
            plot(theta_plot, squeeze(K_all(i,j,:,r)), ...
                'Color',colors(r,:), 'LineWidth',1.8);
        end
        % region boundaries
        for rr = 1:size(Theta1,1)
            xline(Theta1(rr,1),'--k','Alpha',0.35);
            xline(Theta1(rr,2),'--k','Alpha',0.35);
        end
        title(sprintf('K_{%d%d}(\\theta)',i,j),'FontSize',11);
        xlabel('\theta');
        if j==1, ylabel(sprintf('row %d',i)); end
        set(gca,'FontSize',9);
    end
end
sgtitle(sprintf('Switching LPV State-Feedback Gains   (\\gamma_{global} = %.4f)',Gamopt), ...
        'FontSize',13,'FontWeight','bold');
legend(arrayfun(@(r)sprintf('Region %d',r),1:regnum,'Uni',0), ...
       'Location','southoutside','Orientation','horizontal');

%% -------------------- 3. Lyapunov eigenvalues --------------------
figure('Name','Lyapunov Eigenvalues','Color','w','Position',[80 80 1100 480]);
subplot(1,2,1); hold on; grid on; box on;
for r = 1:regnum
    plot(theta_plot, X_eigs(:,:,r)', 'Color',colors(r,:), 'LineWidth',1.6);
end
for rr = 1:size(Theta1,1)
    xline(Theta1(rr,1),'--k','Alpha',0.3);
    xline(Theta1(rr,2),'--k','Alpha',0.3);
end
yline(0,'r-','LineWidth',1.2);
title('Eigenvalues of X(\theta)  (local performance Lyapunov)');
xlabel('\theta'); ylabel('\lambda(X)');
legend(arrayfun(@(r)sprintf('Reg %d',r),1:regnum,'Uni',0),'Location','best');

subplot(1,2,2); hold on; grid on; box on;
for r = 1:regnum
    plot(theta_plot, Xhat_eigs(:,:,r)', 'Color',colors(r,:), 'LineWidth',1.6);
end
for rr = 1:size(Theta1,1)
    xline(Theta1(rr,1),'--k','Alpha',0.3);
    xline(Theta1(rr,2),'--k','Alpha',0.3);
end
yline(0,'r-','LineWidth',1.2);
title('Eigenvalues of Xhat(\theta)  (global stability Lyapunov)');
xlabel('\theta'); ylabel('\lambda(Xhat)');
sgtitle('Parameter-dependent Lyapunov Matrices – All eigenvalues must stay positive','FontSize',12);

%% -------------------- 4. Parameter space partitioning --------------------
figure('Name','Parameter Partitioning','Color','w','Position',[120 120 700 320]);
hold on; box on;
yl = [0.2 0.8];
for r = 1:regnum
    fill([Theta1(r,1) Theta1(r,2) Theta1(r,2) Theta1(r,1)], ...
         [yl(1) yl(1) yl(2) yl(2)], colors(r,:), 'FaceAlpha',0.35, 'EdgeColor','k');
    text(mean(Theta1(r,:)), mean(yl), sprintf('Region %d\n\\gamma=%.3f',r,gamVec(r)), ...
         'HorizontalAlignment','center','FontWeight','bold');
end
xlim([thetaMin-0.05 thetaMax+0.05]);
ylim([0 1]);
set(gca,'YTick',[]);
xlabel('Scheduling parameter \theta');
title(sprintf('Partitioning of parameter space  (\\delta = %g)',delta));
grid on;

%% -------------------- 5. Closed-loop poles at several test points --------------------
test_theta = linspace(thetaMin, thetaMax, 9);
fprintf('\n----- Closed-loop poles at sample points -----\n');
figure('Name','Closed-loop Poles','Color','w','Position',[100 100 900 500]);
hold on; grid on; box on; axis equal;
theta_for_poles = [];
poles_all = [];

for k = 1:length(test_theta)
    th = test_theta(k);
    % find which region
    r = find(th >= Theta1(:,1) & th <= Theta1(:,2), 1, 'first');
    if isempty(r), r = regnum; end
    
    Ft  = F_theta(th);
    Lth = zeros(nu,n);
    for j = 1:nF
        Lth = Lth + Ft(j)*Kopt.L(:,:,j,r);
    end
    Kth = Lth / Kopt.G(:,:,r);          % K = L * inv(G)
    
    % plant at this theta
    Ga  = AugPltEval(Gasym, th);
    Acl = Ga.A + Ga.B2 * Kth;
    p   = eig(Acl);
    
    plot(real(p), imag(p), 'o', 'MarkerSize',8, ...
         'MarkerFaceColor',colors(r,:), 'MarkerEdgeColor','k');
    fprintf('θ = %6.3f (Reg %d): poles = ', th, r);
    fprintf('%.3f%+.3fi  ', [real(p) imag(p)].');
    fprintf('\n');
    
    theta_for_poles = [theta_for_poles; th*ones(n,1)];
    poles_all = [poles_all; p];
end
xline(0,'r-','LineWidth',1.2);
yline(0,'k-');
xlabel('Real part'); ylabel('Imaginary part');
title('Closed-loop poles at 9 sample points (color = region)');
legend(arrayfun(@(r)sprintf('Region %d',r),1:regnum,'Uni',0),'Location','best');

%% -------------------- 6. Time-domain simulation --------------------
% We simulate the closed-loop discrete-time system for three different
% scheduling trajectories to see the behaviour of the states.

Nsim = 80;                          % simulation steps
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
sgtitle('Closed-loop free response under three different scheduling trajectories','FontSize',13);

%% -------------------- 7. Numerical summary of controllers --------------------
fprintf('\n----- Controller data (G and L for each region) -----\n');
for r = 1:regnum
    fprintf('\n========== Region %d (theta in [%.2f , %.2f]) ==========\n', ...
            r, Theta1(r,1), Theta1(r,2));
    fprintf('G =\n'); disp(Kopt.G(:,:,r));
    for j = 1:nF
        fprintf('L(:,:,%d) =\n',j); disp(Kopt.L(:,:,j,r));
    end
    % sample gain at the middle of the region
    th_mid = mean(Theta1(r,:));
    Ft = F_theta(th_mid);
    Lth = zeros(nu,n);
    for j = 1:nF
        Lth = Lth + Ft(j)*Kopt.L(:,:,j,r);
    end
    Kmid = Lth / Kopt.G(:,:,r);
    fprintf('K(theta=%.3f) =\n', th_mid); disp(Kmid);
end

%% -------------------- 8. Quick positivity check --------------------
fprintf('\n----- Quick positivity check of Lyapunov matrices -----\n');
min_eig_X    = min(X_eigs(:));
min_eig_Xhat = min(Xhat_eigs(:));
fprintf('Smallest eigenvalue of X(theta)    over the grid : %.4e\n', min_eig_X);
fprintf('Smallest eigenvalue of Xhat(theta) over the grid : %.4e\n', min_eig_Xhat);
if min_eig_X > 0 && min_eig_Xhat > 0
    fprintf('>>> Both Lyapunov matrices stay positive definite on the grid.\n');
else
    fprintf('>>> WARNING: some eigenvalues are non-positive!\n');
end

fprintf('\n================================================================\n');
fprintf('Visualization finished.  Please send me the figures / key numbers.\n');
fprintf('================================================================\n');


%% Check Lyapunov positivity ONLY at the LMI grid points
fprintf('\n----- Eigenvalues at actual LMI points -----\n');
for r = 1:regnum
    regid1 = REGID(r,2);
    thetaT = ThetaT{regid1};
    for id = 1:length(thetaT)
        th = thetaT(id);
        Ft = F_theta(th);
        Xth  = zeros(n);
        Xhth = zeros(n);
        for j = 1:nF
            Xth  = Xth  + Ft(j)*Kopt.X(:,:,j,r);
            Xhth = Xhth + Ft(j)*Kopt.Xhat(:,:,j,r);
        end
        evX  = sort(real(eig(Xth)),'descend');
        evXh = sort(real(eig(Xhth)),'descend');
        fprintf('Reg %d, θ = %7.3f :  min(eig(X)) = %10.3e ,  min(eig(Xhat)) = %10.3e\n', ...
                r, th, min(evX), min(evXh));
    end
end