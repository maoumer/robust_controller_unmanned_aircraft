% MCEN 6228 Final Project - Mohammed Adib Oumer
%
% The following rigid body model of a rotary-wing MAV has 11 states, 3 inputs, 
% 5 (potential) outputs for controller synthesis. Tracking problem. 
%
% The states (in order) are forward speed (u), lateral velocity (v), roll rate (p),  
% pitch rate (q), roll attitude (phi), pitch attitude (theta), a_s, b_s, heave (w),
% yaw rate (r), yaw rate_feedback (r_fb). The control inputs are lateral and 
% longitudinal cyclic for actuating roll and pitch dynamics (delta_lat, delta_long), 
% and tail rotor (delta_ped) for actuating yaw dynamics. The outputs available 
% for feedback are roll attitude (phi), pitch (theta), and washed out yaw rate 
% (r_fb). You are allowed to consider the pitch rate (q) and roll rate (p) if you 
% are having a difficult time meeting design requirements.

close all; clear all; clc

% Rotorcraft dynamics
A = [-0.1778,zeros(1,4),-9.7807,-9.7807,zeros(1,4);
     0,-0.3104,0,0,9.7807,0,0,9.7807,zeros(1,3);
     -0.3326,-0.5353,zeros(1,4),75.7640,343.86,zeros(1,3);
     0.1903,-0.294,zeros(1,4),172.62,-59.958,zeros(1,3);
     0,0,1,zeros(1,8);zeros(1,3),1,zeros(1,7);
     zeros(1,3),-1,0,0,-8.1222,4.6535,zeros(1,3);
     0,0,-1,zeros(1,3),-0.0921,-8.1222,zeros(1,3);
     zeros(1,6),17.168,7.1018,-0.6821,-0.1070,0;
     0,0,-0.2834,zeros(1,5),-0.1446,-5.5561,-36.674;
     zeros(1,9),2.7492,-11.1120];
 
B = [zeros(6,3);
     0.0632,3.339,0;            % delta_lat (roll motions)
     3.1739,0.2216,0;           % delta_long (pitch motions)
     zeros(1,3);
     0,0,-74.364;               % delta_ped (yaw motions)
     zeros(1,3)];

C = [zeros(1,4),1,zeros(1,6);   % phi
     zeros(1,5),1,zeros(1,5);   % theta
     zeros(1,9),0,1];           % r_fb

D = zeros(size(C,1),size(B,2));

G = ss(A,B,C,D);

%% b)

s = tf('s');
% Performance weighting function
A = 0.005; M = 2; wb = 1;
wp = (s/M+wb)/(s+wb*A); Wp = wp*eye(size(G,1));

% Generalized plant P (plant input tracking configuration)
systemnames = 'G Wp';
inputvar = '[w(3); u(3)]'; % These are the inputs to P
outputvar = '[Wp; w-G]'; % These are the outputs from P
input_to_G = '[u]';
input_to_Wp = '[w-G]';
sysoutname = 'P';
sysic;
P = minreal(ss(P));

% Compute H-infinity controller
n_meas = 3; n_ctrl = 3;
[K,CL,GAM,info] = hinfsyn(P,n_meas,n_ctrl,'method','ric','Tolgam',1e-3,'DISPLAY','on');

% Plots of open loop TF singular values
figure(); clf;
sigma(G,'b',G*K,'r--'); grid on; %,{1e-3,1e2}
legend('\sigma(G)','\sigma(G*K)','Location','best');
title('Open loop singular values');

% Plots of closed loop TF singular values
So = inv(eye(3)+G*K);
To = G*K*inv(eye(3)+G*K);
figure(); clf;
sigma(So,'b',To,'g',inv(Wp),'r--'); grid on; %,{1e-3,1e2}
legend('\sigma(So)','\sigma(To)','\sigma(Wp^{-1})','Location','best');
title('Closed loop singular values');

% Step responses
% The transfer function w->y is To, do->y is So. 
% Row 1 for phi, Row 2 for theta. Column 1 is r_1, column 2 is r_2. r_3 is
% always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(To(1,1:2),100); grid on; %G(1,1:end),'b', - unstable
ylabel('r \rightarrow \phi');
title('Closed loop step response');
subplot(2,1,2), step(To(2,1:2),100); grid on; %G(2,1:end),'b', - unstable
ylabel('r \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(To(1,1),100); grid on;
ylabel('r_1 \rightarrow \phi');
title('Closed loop step response');
subplot(1,2,2), step(To(2,2),100); grid on;
ylabel('r_2 \rightarrow \theta');

% Row 1 for phi, Row 2 for theta. Column 1 is do_1, column 2 is do_2. do_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(So(1,1:2),100); grid on;
ylabel('d_o \rightarrow \phi');
title('Disturbance step response');
subplot(2,1,2), step(So(2,1:2),100); grid on;
ylabel('d_o \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(So(1,1),100); grid on;
ylabel('d_{o_1} \rightarrow \phi');
title('Disturbance step response');
subplot(1,2,2), step(So(2,2),100); grid on;
ylabel('d_{o_2} \rightarrow \theta');

%% COMMENTS
% 
%  Comment on the performance of design in terms of its bandwidth, tracking
%  error and disturbance rejection capabilities.
% 

%% c)

%i) 1st iteration
% Uncertainty weights
wi = (s+0.2)/(0.5*s+1); Wi = wi*eye(3);
% Performance weighting function
A = 0.005; M = 2; wb = 1;
wp = (s/M+wb)/(s+wb*A); Wp = wp*eye(3);

% % Generalized plant P (plant input tracking + uncertainty configuration)
% % Multiplicative input uncertainty 
% systemnames = 'G Wp Wi';
% inputvar = '[ydel(3); w(3); u(3)]'; % These are the inputs to P
% outputvar = '[Wi; Wp; w-G]'; % These are the outputs from P
% input_to_G = '[u+ydel]';
% input_to_Wp = '[w-G]';
% input_to_Wi = '[u]';
% sysoutname = 'P';
% sysic;
% P = minreal(ss(P));
% 
% % Compute H-infinity controller
% n_meas = 3; n_ctrl = 3;
% [K,CL,GAM,info] = hinfsyn(P,n_meas,n_ctrl,'method','ric','Tolgam',1e-3,'DISPLAY','on');

So = inv(eye(3)+G*K);
To = G*K*inv(eye(3)+G*K);
Ti = K*G*inv(eye(3)+K*G);
N = [-Wi*Ti Wi*K*So; -Wp*So*G Wp*So]; N11 = Wi*Ti;

% Plot options
LinMagopt = bodeoptions;
LinMagopt.PhaseVisible = 'off';
LinMagopt.MagUnits = 'abs';
LinMagopt.XLim = [1e-2,1e3];
omega = logspace(-2,3);

% Compute and plot Mu bounds, RS
WiTi_frd = frd(N11,omega);
clear BlkStruct; BlkStruct = [1 0;1 0;1 0]; %diagonal structured complex
[mu_WiTi] = mussv(WiTi_frd, BlkStruct);
figure(); clf;
bodeplot(mu_WiTi(1,1),'*',mu_WiTi(1,2),'r-',LinMagopt);
title('RS \mu plot - diagonal complex'); grid on;

% Compute and plot Mu bounds, RP
N_frd = frd(N,omega);
clear BlkStruct; BlkStruct = [1 0;1 0;1 0;3 3]; %diagonal structured complex
[mu_N] = mussv(N_frd, BlkStruct);
figure(); clf;
bodeplot(mu_N(1,1),'*',mu_N(1,2),'r-',LinMagopt);
title('RP \mu plot - diagonal complex'); grid on;

% % Step responses w->phi, w->theta
% % The transfer function w->y is To, do->y is So.
% % Row 1 for phi, Row 2 for theta. Column 1 is r_1, column 2 is r_2. r_3
% % is always 0 so I didn't include it
% figure(); clf;
% subplot(2,1,1), step(To(1,1:2),100); grid on; %G(1,1:end),'b', - unstable
% ylabel('r \rightarrow \phi');
% title('Closed loop step response');
% subplot(2,1,2), step(To(2,1:2),100); grid on; %G(2,1:end),'b', - unstable
% ylabel('r \rightarrow \theta');
% 
% % Row 1 for phi, Row 2 for theta. Column 1 is do_1, column 2 is do_2. do_3
% % is always 0 so I didn't include it
% figure(); clf;
% subplot(2,1,1), step(So(1,1:2),100); grid on;
% ylabel('d_o \rightarrow \phi');
% title('Disturbance step response');
% subplot(2,1,2), step(So(2,1:2),100); grid on;
% ylabel('d_o \rightarrow \theta');

%% COMMENTS
% 
%  Calculate margins 1/mu_peak on RS and RP. mu_peak_RS = 1.05, mu_peak_RP
%  = 1.55
% 

%%

%ii) 2nd iteration with new Wi, Wp
% Uncertainty weights
wi = (s+0.2)/(0.6*s+1); Wi = eye(3)*wi;
% Performance weighting function
A = 0.005; M = 6; wb = 0.3;
wp = (s/M+wb)/(s+wb*A); Wp = wp*eye(3);

% % Generalized plant P (plant input tracking + uncertainty configuration)
% systemnames = 'G Wp Wi';
% inputvar = '[ydel(3); w(3); u(3)]'; % These are the inputs to P
% outputvar = '[Wi; Wp; w-G]'; % These are the outputs from P
% input_to_G = '[u+ydel]';
% input_to_Wp = '[w-G]';
% input_to_Wi = '[u]';
% sysoutname = 'P';
% sysic;
% P = minreal(ss(P));
% 
% % Compute H-infinity controller
% n_meas = 3; n_ctrl = 3;
% [K,CL,GAM,info] = hinfsyn(P,n_meas,n_ctrl,'method','ric','Tolgam',1e-3,'DISPLAY','on');

So = inv(eye(3)+G*K);
To = G*K*inv(eye(3)+G*K);
Ti = K*G*inv(eye(3)+K*G);
N = [-Wi*Ti Wi*K*So; -Wp*So*G Wp*So]; N11 = Wi*Ti;

% Plot options
LinMagopt = bodeoptions;
LinMagopt.PhaseVisible = 'off';
LinMagopt.MagUnits = 'abs';
LinMagopt.XLim = [1e-2,1e3];
omega = logspace(-2,3);

% Compute and plot Mu bounds, RS
WiTi_frd = frd(N11,omega);
clear BlkStruct; BlkStruct = [1 0;1 0;1 0]; %diagonal structured complex
[mu_WiTi] = mussv(WiTi_frd, BlkStruct);
figure(); clf;
bodeplot(mu_WiTi(1,1),'*',mu_WiTi(1,2),'r-',LinMagopt);
title('RS \mu plot - diagonal complex'); grid on;

% Compute and plot Mu bounds, RP
N_frd = frd(N,omega);
clear BlkStruct; BlkStruct = [1 0;1 0;1 0;3 3]; %diagonal structured complex
[mu_N] = mussv(N_frd, BlkStruct);
figure(); clf;
bodeplot(mu_N(1,1),'*',mu_N(1,2),'r-',LinMagopt);
title('RP \mu plot - diagonal complex'); grid on;

% Step responses w->phi, w->theta
% The transfer function w->y is To, do->y is So.
% Row 1 for phi, Row 2 for theta. Column 1 is r_1, column 2 is r_2. r_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(To(1,1:2),100); grid on; %G(1,1:end),'b', - unstable
ylabel('r \rightarrow \phi');
title('Closed loop step response');
subplot(2,1,2), step(To(2,1:2),100); grid on; %G(2,1:end),'b', - unstable
ylabel('r \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(To(1,1),100); grid on;
ylabel('r_1 \rightarrow \phi');
title('Closed loop step response');
subplot(1,2,2), step(To(2,2),100); grid on;
ylabel('r_2 \rightarrow \theta');

% Row 1 for phi, Row 2 for theta. Column 1 is do_1, column 2 is do_2. do_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(So(1,1:2),100); grid on;
ylabel('d_o \rightarrow \phi');
title('Disturbance step response');
subplot(2,1,2), step(So(2,1:2),100); grid on;
ylabel('d_o \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(So(1,1),100); grid on;
ylabel('d_{o_1} \rightarrow \phi');
title('Disturbance step response');
subplot(1,2,2), step(So(2,2),100); grid on;
ylabel('d_{o_2} \rightarrow \theta');

%% COMMENTS
% 
%  Calculate and comment on improved margins 1/mu_peak on RS and RP. Assess
%  how much uncertainty your design can tolerate. Comment on the
%  performance of design in terms of its bandwidth, tracking error and
%  disturbance rejection capabilities. mu_peak_RS = 0.957, mu_peak_RP =
%  1.12
% 

%% d)

%Uncertainty weights
wi = (s+0.2)/(0.5*s+1); Wi = eye(3)*wi;
% Performance weighting function
A = 0.005; M = 6; wb = 0.3; %6, 0.3
wp = (s/M+wb)/(s+wb*A); Wp = wp*eye(3);

% Generalized plant P (plant input tracking + uncertainty configuration)
% Multiplicative input uncertainty 
systemnames = 'G Wp Wi';
inputvar = '[ydel(3); w(3); u(3)]'; % These are the inputs to P
outputvar = '[Wi; Wp; w-G]'; % These are the outputs from P
input_to_G = '[u+ydel]';
input_to_Wp = '[w-G]';
input_to_Wi = '[u]';
sysoutname = 'P';
sysic;
P = minreal(ss(P));

%%

%Automatic D-K iteration
Delta = [ultidyn('D_1',[1 1]) 0 0; 0 ultidyn('D_2',[1 1]) 0; 0 0 ultidyn('D_3',[1 1])];
Punc = lft(Delta,P); % Close upper LFT
[Kauto,clp,dkinfo] = musyn(Punc,n_meas,n_ctrl);

% D-step
Nf = frd(lft(P,Kauto),omega);
[mu_Nf,mu_info] = mussv(Nf,BlkStruct);
mu_RP = norm(mu_Nf(1,1),inf,1e-6)

% Plot mu
figure(); clf;
bodeplot(mu_Nf(1,1),LinMagopt); grid on;
title('Automatic D-K iteration - musyn');

% Performance plots for final controller
% Closed loop TF
So = inv(eye(3)+G*Kauto);
To = G*Kauto*inv(eye(3)+G*Kauto);

% Open loop singular values
figure(); clf;
sigma(G,'b',G*Kauto,'r--'); grid on;
legend('\sigma(G)','\sigma(G*K)','Location','best');
title('Open Loop Singular Values - Auto');

% Closed loop singular values and wp
figure(); clf;
sigma(So,To,'g-',inv(Wp),'r--'); grid on;
legend('\sigma(S_o)','\sigma(T_o)','\sigma(Wp^{-1})','Location','best');
title('Closed Loop Singular Values - Auto');

% Step responses w->phi, w->theta
% The transfer function w->y is To, do->y is So.
% Row 1 for phi, Row 2 for theta. Column 1 is r_1, column 2 is r_2. r_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(To(1,1:2)); grid on; %G(1,1:end),'b', - unstable
ylabel('r \rightarrow \phi');
title('Closed loop step response - Auto');
subplot(2,1,2), step(To(2,1:2)); grid on; %G(2,1:end),'b', - unstable
ylabel('r \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(To(1,1),35); grid on;
ylabel('r_1 \rightarrow \phi');
title('Closed loop step response');
subplot(1,2,2), step(To(2,2),35); grid on;
ylabel('r_2 \rightarrow \theta');

% Row 1 for phi, Row 2 for theta. Column 1 is do_1, column 2 is do_2. do_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(So(1,1:2)); grid on;
ylabel('d_o \rightarrow \phi');
title('Disturbance step response - Auto');
subplot(2,1,2), step(So(2,1:2)); grid on;
ylabel('d_o \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(So(1,1),35); grid on;
ylabel('d_{o_1} \rightarrow \phi');
title('Disturbance step response');
subplot(1,2,2), step(So(2,2),35); grid on;
ylabel('d_{o_2} \rightarrow \theta');

%%

%Manual D-K iteration
% Initialize
omega = logspace(-3,3);
d0 = 1; D = append(d0,d0,d0,tf(eye(2)),tf(eye(2)),tf(eye(2)));

% Initial K-step: Compute initial H-infinity controller
n_meas = 3; n_ctrl = 3;
[K,CL,GAM(1),info] = hinfsyn(D*P*inv(D),n_meas,n_ctrl,'method','ric','Tolgam',1e-5,'DISPLAY','on');

% Initial D-step: Compute RP level
Nf = frd(lft(P,K),omega);
[mu_Nf,mu_info] = mussv(Nf,BlkStruct);
mu_RP(1) = norm(mu_Nf(1,1),inf,1e-6)

% Plot options
LinMagopt = bodeoptions;
LinMagopt.PhaseVisible = 'off';
LinMagopt.MagUnits = 'abs';
LinMagopt.XLim = [1e-2,1e3];
omega = logspace(-2,3);

% Generate mu plot
figure(); clf;
bodeplot(mu_Nf(1,1),LinMagopt);
title('Manual D-K iteration'); grid on; hold on;

% D-K iteration loop
n = 6; %past this iteration, I get controllability/observability error for D*P*inv(D)
for i = 2:n
    
    %Fit resulting D-scales
    [dsysl,dsysr] = mussvunwrap(mu_info);
    dsysl = dsysl/dsysl(3,3);
    di = fitfrd(genphase(dsysl(1,1)),3);
    Di = append(di,di,di,tf(eye(2)),tf(eye(2)),tf(eye(2)));
    
    % K-step: Compute H-infinity controllers
    [Ki,CL,GAM(i),info] = hinfsyn(Di*P*inv(Di),n_meas,n_ctrl,'method','ric','Tolgam',1e-3,'DISPLAY','on');
    
    % D-step: Compute RP level
    Nf = frd(lft(P,Ki),omega);
    [mu_Nf,mu_info] = mussv(Nf,BlkStruct);
    mu_RP(i) = norm(mu_Nf(1,1),inf,1e-6)
    
    % Add to plot
    if i == n
        bodeplot(mu_Nf(1,1),LinMagopt); grid on;
        legend('First iteration','Final iteration');
        Kfinal = Ki;
    end
end

% Peak of mu plot as a function of iteration
figure(); clf;
plot(mu_RP,'b-*'); hold on;
plot(GAM,'r--*'); grid on;
xlabel('Iteration');
ylabel('\mu and \gamma');
title('RP level and H_{\infty} cost');
legend('Peak of \mu','\gamma');

% Performance plots for final controller
% Closed loop TF
So = inv(eye(3)+G*Kfinal);
To = G*Kfinal*inv(eye(3)+G*Kfinal);

% Open loop singular values
figure(); clf;
sigma(G,'b',G*Kfinal,'r--'); grid on;
legend('\sigma(G)','\sigma(G*K)','Location','best');
title('Open Loop Singular Values - Manual');

% Closed loop singular values and wp
figure(); clf;
sigma(So,To,'g-',inv(Wp),'r--'); grid on;
legend('\sigma(S_o)','\sigma(T_o)','\sigma(Wp^{-1})','Location','best');
title('Closed Loop Singular Values - Manual');

% Step responses w->phi, w->theta
% The transfer function w->y is To, do->y is So.
% Row 1 for phi, Row 2 for theta. Column 1 is r_1, column 2 is r_2. r_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(To(1,1:2)); grid on; %G(1,1:end),'b', - unstable
ylabel('r \rightarrow \phi');
title('Closed loop step response - Manual');
subplot(2,1,2), step(To(2,1:2)); grid on; %G(2,1:end),'b', - unstable
ylabel('r \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(To(1,1),35); grid on;
ylabel('r_1 \rightarrow \phi');
title('Closed loop step response');
subplot(1,2,2), step(To(2,2),35); grid on;
ylabel('r_2 \rightarrow \theta');

% Row 1 for phi, Row 2 for theta. Column 1 is do_1, column 2 is do_2. do_3
% is always 0 so I didn't include it
figure(); clf;
subplot(2,1,1), step(So(1,1:2)); grid on;
ylabel('d_o \rightarrow \phi');
title('Disturbance step response - Manual');
subplot(2,1,2), step(So(2,1:2)); grid on;
ylabel('d_o \rightarrow \theta');

figure(); clf;
subplot(1,2,1), step(So(1,1),35); grid on;
ylabel('d_{o_1} \rightarrow \phi');
title('Disturbance step response');
subplot(1,2,2), step(So(2,2),35); grid on;
ylabel('d_{o_2} \rightarrow \theta');

%% COMMENTS
% 
%  Calculate and comment on improved RP and improved margins 1/mu_peak on
%  RS and RP. Comment on the performance of final design in terms of its
%  bandwidth, tracking error and disturbance rejection capabilities, in
%  particular tradeoff between high bandwidth and robust performance in the
%  presence of uncertainty. mu_peak_RP = 0.816 (auto),mu_peak_RP = 0.983
%  (manual)
% 
