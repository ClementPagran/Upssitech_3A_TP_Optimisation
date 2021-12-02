clear all
close all

A_REMPLIR = nan;

vR = 1, omegaR = 0.2, Te = 0.25

Ac = [0 omegaR 0 ; -omegaR 0 vR ; 0 0 0]; Bc = [-1 0 ; 0 0 ; 0 -1];
Ad = ...
[  cos(Te*omegaR), sin(Te*omegaR), (vR - vR*cos(Te*omegaR))/omegaR;
  -sin(Te*omegaR), cos(Te*omegaR),     (vR*sin(Te*omegaR))/omegaR;
              0,            0,                          1];
Bd = ...
[          -sin(Te*omegaR)/omegaR, (vR*(sin(Te*omegaR) - Te*omegaR))/omegaR^2;
  (2*sin((Te*omegaR)/2)^2)/omegaR,    -(2*vR*sin((Te*omegaR)/2)^2)/omegaR^2;
                             0,                                   -Te];
Gc = ss(Ac,Bc,eye(3),0); clear Ac Bc
B0Gd = ss(Ad,Bd,eye(3),0,Te); B0Gd_bis = c2d(Gc,Te,'zoh'); 
     
%% Section IV - DLQ Riccati
fprintf('= Section IV - DLQ Riccati =====================================\n');
N = 100; % Attention, horizon de taille N => indices de 1 à N+1 !
nb_cas = 3;
S = {100*eye(3),100*eye(3),diag([1 1 100])};
Q = {100*eye(3),1*eye(3),diag([100 100 1])};
R = {1*eye(2),100*eye(2),diag([1 10])};
TIME = [0:N]*Te;

%
% Résolution de l'équation récurrente rétrograde de Riccati
P = {nan(3,3,N+1),nan(3,3,N+1),nan(3,3,N+1)};
for i_cas = 1:nb_cas
    P{i_cas}(:,:,N+1) = S{i_cas}; % remplissage de P[N]
    for i_time = N+1:-1:2
        % depuis les instants k=N jusqu'à k=1 afin de remplir P[N-1] à P[0]
        P{i_cas}(:,:,i_time-1) = Ad' * P{i_cas}(:,:,i_time) * Ad + Q{i_cas} - Ad' * P{i_cas}(:,:,i_time) * Bd * inv(R{i_cas} + Bd' * P{i_cas}(:,:,i_time) * Bd) * (Bd' * P{i_cas}(:,:,i_time) * Ad);
    end
end
%
% Gain de retour d'état u[k] = -K[k] x[k] depuis k=N-1 jusqu'à k=0
K = {nan(2,3,N+1),nan(2,3,N+1),nan(2,3,N+1)}; % dernière page = NaN's
for i_cas = 1:nb_cas
    for i_time = 1:N
        K{i_cas}(:,:,i_time) = inv(R{i_cas} + Bd' * P{i_cas}(:,:,i_time+1) * Bd) * Bd' * P{i_cas}(:,:,i_time+1) * Ad;
    end
end
%
% Simulation des états depuis k=0 jusqu'à k=N
U = {nan(2,1,N+1),nan(2,1,N+1),nan(2,1,N+1)};
X = {nan(3,1,N+1),nan(3,1,N+1),nan(3,1,N+1)};
for i_cas = 1:nb_cas
    X{i_cas}(:,:,1) = [4 2 pi/3]';
    for i_time = 1:N
        U{i_cas}(:,:,i_time) = -K{i_cas}(:,:,i_time) * X{i_cas}(:,:,i_time);
        X{i_cas}(:,:,i_time+1) = Ad * X{i_cas}(:,:,i_time) +  Bd * U{i_cas}(:,:,i_time);
    end
end
%
% Plot
for i_cas = 1:nb_cas
    plotX{i_cas}=squeeze(X{i_cas});
    plotU{i_cas}=squeeze(U{i_cas});
    figure(i_cas); title(sprintf('== DLQR -- Cas %d ==',i_cas));
    subplot(2,1,1); title('x');
    plot(TIME,plotX{i_cas}(1,:),'x-',...
         TIME,plotX{i_cas}(2,:),'x-.',...
         TIME,plotX{i_cas}(3,:),'x:');
    legend('x1','x2','x3');
    subplot(2,1,2); title('u');
    plot(TIME,plotU{i_cas}(1,:),'x-',...
         TIME,plotU{i_cas}(2,:),'x-.');
    legend('u1','u2');
end

%% Section VI - DLQR
fprintf('= Section VI - DLQR =====================================\n');
for i_cas = 1:nb_cas
    fprintf('Cas %d\n',i_cas);
    fprintf('--> solutions DLQR\n');
    P_DLQR{i_cas} = dlqr(Ad,Bd,Q{i_cas},R{i_cas});
    K_DLQR{i_cas} = inv(R{i_cas} + Bd' * P{i_cas}(:,:,i_time) * Bd) * Bd' * P{i_cas}(:,:,i_time) * Ad;
    P_DLQR{i_cas}, K_DLQR{i_cas}
    fprintf('--> "régimes permanents" P et K du problème DLQ\n');
    A_REMPLIR;
end

%% Section V - DLQ par Programmation Quadratique
fprintf('= Section V - DLQ par Programmation Quadratique =====================================\n');
for i_cas = 1:nb_cas
    barQ{i_cas} = kron(eye(N),Q;
    barR{i_cas} = A_REMPLIR;
    %
    barS = A_REMPLIR;
    barT = A_REMPLIR;
    for i_barST = 1:N
        barS = A_REMPLIR;
        barT = A_REMPLIR;
    end
    H = A_REMPLIR;
    f = A_REMPLIR;
    Y = A_REMPLIR;
    %
    xiU{i_cas} = A_REMPLIR;
    clear barS barT H f Y
    fprintf('Cas %d\n',i_cas);
    fprintf('--> u_B^*[0:N-1] solution DLQ Riccati\n');
    plotU{i_cas}(:,1:(end-1))
    fprintf('--> u_B^*[0:N-1] solution Prog Quad\n');
    reshape(xiU{i_cas},2,N)
    fprintf('--> Écart absolu maxi\n');
    max(max(abs(plotU{i_cas}(:,1:(end-1))-reshape(xiU{i_cas},2,N))))
end

%% Section VII - DLQ avec contraintes par Programmation Quadratique
fprintf('= Section V - DLQ avec contraintes par Programmation Quadratique =====================================\n');
i_cas = 1; % On choisit le Cas 1
%
% Critère 1/2 xi^T H xi où xi = [E^T U^T]^T
% avec E = [e[1];...;e[N]] et U = [uB[1];...;uB[N-1]]
H = A_REMPLIR;
f = A_REMPLIR;
%
% Contraintes égalités (indépendantes du cas considéré)
Ae = A_REMPLIR;
Be = A_REMPLIR;
Ae1 = A_REMPLIR;
Ae2 = A_REMPLIR;
Ae = [Ae1 Ae2];
clear Ae1 Ae2
Be(1:3,:) = A_REMPLIR;
%
% Contraintes inégalités -- on considèrera nb_cas_ineg cas différents
nb_cas_ineg = 1; % on peut aller jusqu'à 3 par exemple
xiEU = cell(nb_cas_ineg);
fval = cell(nb_cas_ineg);
exitflag = cell(nb_cas_ineg);
Ai = cell(nb_cas_ineg);
Bi = cell(nb_cas_ineg);
uMIN = cell(nb_cas_ineg);
uMAX = cell(nb_cas_ineg);
L = cell(nb_cas_ineg); % sorties critiques z = L e
zMIN = cell(nb_cas_ineg);
zMAX = cell(nb_cas_ineg);
uMIN{1} = [-3;-0.5]; uMAX{1} = [3;2];
zMIN{1} = -0.7; zMAX{1} = pi/3; L{1} = [0 0 1]; %contrainte sur e3
for i_cas_ineg = 1:nb_cas_ineg
    Ai{i_cas_ineg} = zeros(2*N+4*N,3*N+2*N);
    Bi{i_cas_ineg} = zeros(2*N+4*N,1);
    Ai{i_cas_ineg} = A_REMPLIR;
    Bi{i_cas_ineg} = A_REMPLIR;
    [xiEU{i_cas_ineg},fval{i_cas_ineg},exitflag{i_cas_ineg}] = ...
        quadprog(H,f,Ai{i_cas_ineg},Bi{i_cas_ineg},Ae,Be);
    if (exitflag{i_cas_ineg} ~= 1), error('Revoir programme optim'); end
    %
    U_in_xiEU = xiEU{i_cas_ineg}((3*N+1):(5*N),:);
    plotUcontraint{i_cas_ineg} = reshape(U_in_xiEU,2,N); clear U_in_xiEU
    plotUcontraint{i_cas_ineg} = [plotUcontraint{i_cas_ineg} nan(2,1)];
    %
    X_in_xiEU = xiEU{i_cas_ineg}(1:(3*N),:);
    plotXcontraint{i_cas_ineg} = reshape(X_in_xiEU,3,N); clear X_in_xiEU
    plotXcontraint{i_cas_ineg} = [X{i_cas}(:,:,1) plotXcontraint{i_cas_ineg}];
end
%
% Plot
for i_cas_ineg = 1:nb_cas_ineg
    figure(10+i_cas_ineg);
    title(sprintf('== DLQR CONTRAINT -- Cas %d -- Cas\\_ineg %d ==',i_cas,i_cas_ineg));
    subplot(2,1,1); title('x');
    plot(TIME,plotXcontraint{i_cas}(1,:),'o-',...
         TIME,plotXcontraint{i_cas}(2,:),'o-.',...
         TIME,plotXcontraint{i_cas}(3,:),'o:');
    legend('x1','x2','x3');
    subplot(2,1,2); title('u');
    plot(TIME,plotUcontraint{i_cas}(1,:),'o-',...
         TIME,plotUcontraint{i_cas}(2,:),'o-.');
    legend('u1','u2');
end
