%% Define mesh (spatial)

eta_l = 0; % lower boundary
eta_t = 4; % top boundary
J = 100; % number of grid points
deta = (eta_t-eta_l)/(J); % vertical step
eta = linspace(eta_l,eta_t,J); % vertical grid points

%% Constants and dimensions

% Dimensional constants
kappa_dim = 0.3; %m^2s^-1
F_dim = 16e-3; %m^2s^-2
N_dim = 2.16e-2; %s^-1
mu_dim = 1e-6; %s^-1
k_dim = 2*pi/(40e6); %m^-1
c_dim = 30; %m/s

H = 7000; %m
Omega = 7.27e-5; % rad/s 
a = 6.371e6; % m 
beta_dim = 2*Omega/a; % s^-1 

chi_dim = 10; %ppmv
L = 1e6; % m
R = 287; % J kg^−1 K^−1 = m^2 s^−2 K^−1
h = @(j) (5.4e-7*(1+(2/3)*j*deta).*(j*deta<=1.99)+1.56e-6.*(j*deta>1.99))/100; % s^-1
gamma_O = -2.24e-7; % s^-1
gamma_T = -1.24e-12*1e6; % m^-1 s
s_O = 8.39e-2*1e-6; % m s^-3
s_C = 0; %-1e-7*1e-6; % m s^-3

% Dimensionless constants
alpha = @(j) (0.55+0.55*(j*deta)).*(j*deta<=1.99)+1.65.*(j*deta>1.99);
Lambda = kappa_dim*N_dim*mu_dim/(k_dim*c_dim*F_dim);
eps = N_dim*mu_dim/(H*k_dim*c_dim^2);
beta = beta_dim/(k_dim^2*c_dim);
theta = 2*pi*k_dim*c_dim^3/(360*N_dim*mu_dim*F_dim);

Gamma_O = gamma_O*k_dim*c_dim^3/(mu_dim*F_dim*N_dim);
Gamma_T = gamma_T*L^2*Omega*c_dim^2/(a*F_dim*chi_dim);

Xi = @(j) h(j)*L^2*Omega*mu_dim/(a*k_dim*F_dim*N_dim);
S_O = s_O*c_dim*chi_dim/(F_dim*N_dim^2);
S_C = s_C*c_dim*chi_dim/(F_dim*N_dim^2);

% Dimensionless variables
zz = @(eta) eta*(k_dim*c_dim^2/(N_dim*mu_dim*1000))+17; %km
tt = @(xi) xi*(k_dim*c_dim^3/(N_dim*mu_dim*F_dim*3600*24)); % days
uu = @(u) c_dim*u; % m/s
TT = @(T) L^2*Omega*H*N_dim*mu_dim/(R*a*k_dim*c_dim)*T; % K
cchi = @(chi) chi_dim*chi; % ppv m^-1

%% Define time step

dt = 0.05; % days time step
dxi = dt*(3600*24)*(N_dim*mu_dim*F_dim)/(k_dim*c_dim^3); % convert
N = 99999; % # timesteps
xi = 0:dxi:dxi*N; % time array

%% Initialise arrays

u = zeros(length(eta),N+1); % mean wind array
T = zeros(length(eta),N+1); % temperature array
chi_O = zeros(length(eta),N+1); % mixing ratio array
w = zeros(length(eta),N+1); % upwelling array

% initial states
u(:,1) = 0.5*sin(pi*eta/4);
T(:,1) = 0.5*(pi/4)*cos(pi*eta/4);
chi_O(:,1) = 0.01*sin(pi*eta/2);
w(:,1) = -Xi(1:J)'.*T(1:J,1)+S_O*chi_O(1:J,1); 

%% Define terms in equation

F_0K = 1; F_0RG = -1; % wave momentum flux
c_K = 1; c_RG = -1; % ms^-1 wavespeeds
k_K = 1; k_RG = 3; % km^-1 wavenumber

g_K = @(u,j,jj,n) alpha(j)./(k_K.*(u(jj,n)-c_K).^2); % rate of decay
g_RG = @(u,j,jj,n) alpha(j)./(k_RG.*(u(jj,n)-c_RG).^2).*(beta./(k_RG^2.*(u(jj,n)-c_RG))-1); % rate of decay

F_K = @(u,j,n) F_0K.*exp(-trapz(g_K(u,j,1:j,n))*deta); % wave momentum flux
dF_K = @(u,j,n) -F_K(u,j,n).*g_K(u,j,j,n); % wave momentum flux derivative

F_RG = @(u,j,n) F_0RG.*exp(-trapz(g_RG(u,j,1:j,n))*deta); % wave momentum flux
dF_RG = @(u,j,n) -F_RG(u,j,n).*g_RG(u,j,j,n); % wave momentum flux derivative

G = @(u,j,n) -exp(eps*deta*j)*(dF_K(u,j,n)+dF_RG(u,j,n)); % forcing term

chi_V = 345/chi_dim;

%% Time stepping

% Create movie
figure
subplot(4,1,1)
plot(uu(u(:,1)),zz(eta),'b--')
hold on
    xline(0,'k:')
hold off
xlim([-30,30])
subplot(4,1,2)
plot(TT(T(:,1)),zz(eta),'r--')
hold on
    xline(0,'k:')
hold off
subplot(4,1,3)
plot(cchi(chi_O(:,1)),zz(eta),'g--')
hold on
    xline(0,'k:')
hold off
subplot(4,1,4)
plot(w(:,1)*F_dim/c_dim,zz(eta),'k--')
hold on
    xline(0,'k:')
hold off
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
M(N) = struct('cdata',[],'colormap',[]);

% Iterate
my_waitbar = waitbar(0,'Calculating QBO... (0%)');
for n = 2:N

    % Update u
    a_0 = @(j) 1 + 2*Lambda*dxi/deta^2 + w(j,n-1)*dxi/deta;
    a_m1 = @(j) - Lambda*dxi/deta^2 - w(j,n-1)*dxi/deta;
    a_p1 = - Lambda*dxi/(deta)^2;
    a_LHS = diag(a_0(2:J-1),0) + diag(a_m1(3:J-1),-1) + diag(a_p1*ones(1,J-3),1);
    
    GG = arrayfun(@(j)G(u,j,(n-1))*(abs(G(u,j,(n-1)))<100)+100*(abs(G(u,j,(n-1)))>=100),2:J-1)';
    a_RHS = u(2:J-1,n-1) + dxi*GG;

    a_LHS = cat(2,[a_m1(2) zeros(1,J-3)]',a_LHS,[zeros(1,J-3) a_p1]');
    a_LHS = cat(1,[1 zeros(1,J-1)],a_LHS,[zeros(1,J-3) 1 -4 3]); % BC LHS
    a_RHS = cat(1,0,a_RHS,0); % BC RHS

    u(1:J,n) = a_LHS\a_RHS; 

    % Update T
    T(1,n) = - (-3*u(1,n)+4*u(2,n)-u(3,n))/(2*deta);
    T(2:J-1,n) = - (u(3:J,n)-u(1:J-2,n))/(2*deta);
    T(J,n) = - (u(J-2,n)-4*u(J-1,n)+3*u(J,n))/(2*deta);

    % Update chi
    b_0 = @(j) 1 + w(j,n-1).*dxi/deta.*(w(j,n-1)>=0) -  w(j,n-1).*dxi/deta.*(w(j,n-1)<0) - dxi*Gamma_O;
    b_m1 = @(j) - w(j,n-1).*dxi/deta.*(w(j,n-1)>=0); % if w>0 backwards difference
    b_p1 = @(j) w(j,n-1).*dxi/deta.*(w(j,n-1)<0); % if w<0 forwards difference

    b_LHS = diag(b_0(2:J-1),0)+diag(b_m1(3:J-1),-1)+diag(b_p1(2:J-2),1);
    b_LHS = cat(2,[b_m1(2) zeros(1,J-3)]',b_LHS,[zeros(1,J-3) b_p1(J-1)]');
    b_RHS = chi_O(2:J-1,n-1)+dxi*Gamma_T*T(2:J-1,n);

    if w(1,n-1) > 0 % BC needed at the beginning
       b_LHS = cat(1,[1 zeros(1,J-1)],b_LHS);
       b_RHS = cat(1,0,b_RHS);
    else % Point calculated using explicit forward differences
       b_LHS = cat(1,[b_0(1) b_p1(1) zeros(1,J-2)],b_LHS);
       b_RHS = cat(1,chi_O(J,n-1)+dxi*Gamma_T*T(J,n),b_RHS);
    end

    if w(J,n-1) < 0 % BC needed at the end
        b_LHS = cat(1,b_LHS,[zeros(1,J-3) 1 -4 3]);
        b_RHS = cat(1,b_RHS,0);
    else % Point calculated using implicit backwards differences
        b_LHS = cat(1,b_LHS,[zeros(1,J-2) b_m1(J-1) b_0(J)]);
        b_RHS = cat(1,b_RHS,chi_O(J,n-1)+dxi*Gamma_T*T(J,n));
    end

    chi_O(1:J,n) = b_LHS\b_RHS;

    % Update w
    w(1:J,n) = -Xi(1:J)'.*T(1:J,n)+S_O*chi_O(1:J,n)+S_C*chi_V*ones(J,1);

    % Update movie
    subplot(4,1,1)
    title("Time: "+string(tt(n*dxi)))
    plot(uu(u(:,n)),zz(eta),'b--')
    hold on
        xline(0,'k:')
    hold off

    subplot(4,1,2)
    plot(TT(T(:,n)),zz(eta),'r--')
    hold on
        xline(0,'k:')
    hold off

    subplot(4,1,3)
    plot(cchi(chi_O(:,n)),zz(eta),'g--')
    hold on
        xline(0,'k:')
    hold off

    subplot(4,1,4)
    plot(w(:,n)*F_dim/c_dim,zz(eta),'k--')
    hold on
        xline(0,'k:')
    hold off

    axis tight manual
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    M(N) = struct('cdata',[],'colormap',[]);
    drawnow
    M(n) = getframe;
    
    % Update waitbar
    waitbar(n/(N-2),my_waitbar,sprintf('Calculating QBO... (Percentage done: %3.3f )',n/(N-2)*100))
end

%% Find period at height 24km

eta_eval = int64((24-17)*(N_dim*mu_dim*1000)/(k_dim*c_dim^2)); % eta corresponding to 24km
s = u(eta_eval/deta,:); % zonal mean flow at 24km
[ac,lags] = xcorr(s,s); % autocorrelation function
locs = islocalmax(ac); % find peaks
period = tt(mean(diff(lags(locs)*dxi)))/30.4167; % period in months

%% Save averages

fileID = fopen('results.txt','w');
fprintf(fileID,'%12s %12s %12s\n','var','mean','std');
fprintf(fileID,'%12s %12f %12f\n','period', period,0);
fprintf(fileID,'%12s %12f %12f\n','u', mean(uu(u(:))),std(uu(u(:))));
fprintf(fileID,'%12s %12f %12f\n','T', mean(TT(T(:))),std(TT(T(:))));
fprintf(fileID,'%12s %12f %12f\n','chi', mean(cchi(chi_O(:))),std(cchi(chi_O(:))));
fprintf(fileID,'%12s %12f %12f\n','w', mean(w(:)*F_dim/c_dim),std(w(:)*F_dim/c_dim));
fclose(fileID);
