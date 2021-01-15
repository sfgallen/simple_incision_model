clear; clc

% model run time
run_time = 7e5;
stabil = 1;
n_plots = 6;

% define and decalte the domain
l_crit = 1e3;
l_max = 1e5;
dx = 100;
L = l_crit:dx:l_max;

% define hacks constants and DA
ka = 6.69;
h = 1.67;
A = 6.69.*L.^1.67;

% define the steam power parameters
Ui = 1e-3;
Ki = 1e-5;

Uf = 3e-3;
Kf = 1e-5;

m = 0.5;
n = 1;

mn = m/n;

%%
% calculate the steady-state slope and elevation from Flint's Law
[S,Lo,Z] = SS_profile(Ui,Ki,m,n,A,L);

% calculate chi
Achi = A.^-mn;
Achi = fliplr(Achi);
chi = zeros(1,length(L));

for i = 2:length(chi)
    dchi = ((Achi(i)+Achi(i-1))/2)*(Lo(i)-Lo(i-1));
    chi(i) = chi(i-1) + dchi;
end
chi = fliplr(chi);

% CFL criterion
% calculate final steady state
[Sf,~,Zf] = SS_profile(Uf,Kf,m,n,A,L);

% calculate dt
dt = min([dx./(Ki.*A.^m.*S.^(n-1)),dx./(Kf.*A.^m.*Sf.^(n-1))])/stabil;

% calculate erosion rate
E = Ki.*A.^m.*S.^n;

% plot initial conditions
subplot(2,2,1)
plot(L,Z,'k-'); hold on
xlabel('Distance (m)'); ylabel('Elevation (m)');
subplot(2,2,2)
plot(A,S,'k-'); hold on
xlabel('Drainage Area'); ylabel('Slope');
set(gca,'xscale','log','yscale','log');
subplot(2,2,3);
plot(L,E,'k-'); hold on
xlabel('Distance (m)'); ylabel('Erosion rate (m/yr)');
subplot(2,2,4);
plot(chi,Z,'k-'); hold on
xlabel('\chi'); ylabel('Elevation');
set(gca,'xdir','reverse');

% calculate erosion rate
Ef = Kf.*A.^m.*Sf.^n;

% plot final conditions
subplot(2,2,1)
plot(L,Zf,'r-'); hold on
xlabel('Distance (m)'); ylabel('Elevation (m)');
subplot(2,2,2)
plot(A,Sf,'r-'); hold on
xlabel('Drainage Area'); ylabel('Slope');
set(gca,'xscale','log','yscale','log');
subplot(2,2,3);
plot(L,Ef,'r-'); hold on
xlabel('Distance (m)'); ylabel('Erosion rate (m/yr)');
subplot(2,2,4);
plot(chi,Zf,'r-'); hold on
xlabel('\chi'); ylabel('Elevation');
set(gca,'xdir','reverse');

% get things step up for the time loop
t_steps = ceil(run_time/dt);
t_plots = floor(t_steps/n_plots);

E_mean = nan(t_steps+1,1);
mod_time = zeros(t_steps+1,1);

E_mean(1) = mean(E);

% set up waitbar
h = waitbar(0,'The model is running...');

% initiate for loop
for t = 1:t_steps
    
    % calculate erosion rate
    E = Kf.*A.^m.*S.^n;
    E_mean(t+1) = mean(E);
    mod_time(t+1) = t*dt;
    
    % evolve the profile
    Z_cur = Z + (Uf.*dt) - (E.*dt);
   
    % set the lower boundary condition
    Z_cur(end) = 0;
    
    % update slope vector
    S = calc_slope(Z_cur,L);
    
    % set Z to Z_cur for next iteration
    Z = Z_cur;
    
    % plot data as needed
    if rem(t,t_plots) == 0
        subplot(2,2,1)
        plot(L,Z,'b-'); hold on
        xlabel('Distance (m)'); ylabel('Elevation (m)');
        subplot(2,2,2)
        plot(A,S,'b-'); hold on
        xlabel('Drainage Area'); ylabel('Slope');
        set(gca,'xscale','log','yscale','log');
        subplot(2,2,3);
        plot(L,E,'b-'); hold on
        xlabel('Distance (m)'); ylabel('Erosion rate (m/yr)');
        subplot(2,2,4);
        plot(chi,Z,'b-'); hold on
        xlabel('\chi'); ylabel('Elevation');
        set(gca,'xdir','reverse');
    end
    waitbar((t/t_steps),h);
end
close(h)

% plot the erosion rate history
figure()
plot(mod_time,E_mean*1000,'b-');
xlabel('Time (yrs)'); ylabel('Erosion rate (mm/yr)');