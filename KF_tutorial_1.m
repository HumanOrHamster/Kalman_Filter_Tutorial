clear;clc;close all

% 1D Position Estimate KF example
% v1 by Binh Bui, 5/16/2017
% Copyright (C) 2017 Binh Bui
% This program is for educational purpose only. 

rng(100)

%% run time
% duration of simulation (sec)
duration = 400;

% time step
dt = .01;
% time vector
t = 0:dt:duration;
n = length(t);

% measurement update interval
updateTime = 1;

%% Linear discrete model setup
A = [1 dt .5*dt^2; 
    0 1 dt; 
    0 0 1];
B = [.5*dt^2 dt 0]';
C = [ 1 0 -.5*dt^2];

%% set up noise parameters

% Noise specification
sigmaMeasureNoise = 20; % position measurement noise (meter)
sigmaAccelNoise = .1; % acceleration noise (m/sec^2)

% Noise Covariance Matrix
Q = B*sigmaAccelNoise^2*B';
% Measurement Error Matrix
R = sigmaMeasureNoise^2; 

% Initial Covariance Matrix
P = diag( [5 1 .1] ).^2

%% Problem initialization

% True States
xTrue = zeros(3,n);
xTrue(:,1) = [0;0;0];
yTrue = zeros(1,n);
yTrue(1) = C*xTrue(:,1);

% Estimated states
xHat = zeros(3,n);
% Initial estimate
xHat(:,1) = [6;1;0.07];
yHat = zeros(1,n);

% Uncorrected states
xUnbound = zeros(3,n);
xUnbound (:,1) = xHat(:,1);

% Accelerometer input vector;
omega = (2*pi)/(duration/5); % frequency of input (arbitrary)
uTrue = 1*sin(omega.*t)*omega; % clean input signal (arbitrary)
constantBias = .1; % bias of the accelerometer 
% add random noise and bias to the input
uNoisy = uTrue+sigmaAccelNoise*randn(1,n)+constantBias; % Noisy input

% Measurement Noise vector
MeasurementNoise = sigmaMeasureNoise*randn(1,n); % measurement noise
yMeasure = zeros(1,n);

% Predefined arrays to store measurements
yFix = []; yFix = [yFix 0];
yFixErr = [];yFixErr = [yFixErr MeasurementNoise(1)];
tFix = []; tFix = [tFix 0];

% std deviation vector initialization
sigma.pos = zeros(1,n);
sigma.pos(1) = sqrt(P(1,1));
sigma.vel = zeros(1,n);
sigma.vel(1) = sqrt(P(2,2));
sigma.bias = zeros(1,n);
sigma.bias(1) = sqrt(P(3,3));
% Error vectors vector initialization
error.pos = zeros(1,n);
error.pos(1) = xTrue(1) - xHat(1);
error.vel = zeros(1,n);
error.vel(1) = xTrue(2) - xHat(2);
error.bias = zeros(1,n);
error.bias(1) = xTrue(3) - xHat(3);

% gain vector initialization
gainStore = zeros(3,n); % store gain value

%% Run KF Simulation

for k =2:n
    % Truth model
    % Truth state propagation
    xTrue(:,k) = A*xTrue(:,k-1)+B*uTrue(k-1);
    % Truth measurement
    yTrue(:,k) = C*xTrue(:,k);
    
    % Estimated model
    % i) Estimated state propagation
    xHat(:,k) = A*xHat(:,k-1)+B*uNoisy(k-1);
    % ii) Estimate measurement
    yHat(k) = C*xHat(:,k);
    
    % Unbounded Propagation (for comparison purpose only)
    xUnbound(:,k) = A*xUnbound(:,k-1)+B*uNoisy(k-1);
    
    % iii) Covariance Propagation
    P = A*P*A' + Q;
    
    % Measurement + Noise
    yMeasure(k) = yTrue(k) + MeasurementNoise(k);
    
    % IF fix is available
    if mod((k-1)*dt,updateTime)==0 
        
        % Store fix data for later plotting
        tFix = [tFix t(k)];
        yFix = [yFix yMeasure(k)];
        yFixErr = [yFixErr MeasurementNoise(k)];
        
        % iv) Residual
        residual = yMeasure(k)-yHat(k);
        
        % v) Gain calculation
        K = P*C'*inv(C*P*C'+R);
        
        % vi) Covariance Update
        P =(eye(3)-K*C)*P;
        
        gainStore(:,k) = K;
        
        
        % vii) State estimate update
        xHat(:,k) = xHat(:,k)+K*residual;
        
    else 
        % no update to xhat
        gainStore(:,k) = gainStore(:,k-1);
    end

    % extract std dev from the diagonals of the Covariance matrix
    sigma.pos(k) = sqrt(P(1,1));
    sigma.vel(k) = sqrt(P(2,2));
    sigma.bias(k) = sqrt(P(3,3));
end
% calulate errors
error.pos = xTrue(1,:) - xHat(1,:);
error.vel = xTrue(2,:) - xHat(2,:);
error.bias = -constantBias.*ones(size(t)) - xHat(3,:);
%%
close all
figure(1)
subplot(3,1,1)
plot(tFix,yFix,'+-','color',[.1 .1 .1],'markersize',5); hold on
plot(t,xTrue(1,:),'b-.','linewidth',1.5); hold on
plot(t,xHat(1,:),'r','linewidth',1.5); hold on
title('Position')
xlabel('sec');ylabel('meter');
% plot(t,xUnbound(1,:),'color',[0 .7 0]) % uncompensated
legend('Measured','True','Estimate')

subplot(3,1,2)
plot(t,xTrue(2,:),'b-.','linewidth',1.5); hold on
plot(t,xHat(2,:),'r','linewidth',1.5); hold on
title('Velocity')
xlabel('sec');ylabel('meter/sec')
% plot(t,xUnbound(2,:),'color',[0 .7 0]) % uncompensated
legend('True','Estimate')

subplot(3,1,3)
plot(t,-constantBias.*ones(size(t)),'b-.','linewidth',1.5); hold on
plot(t,xHat(3,:),'r','linewidth',1.5); hold on
plot(t,xUnbound(3,:),'color',[0 .7 0])
title('Accleration Bias')
xlabel('sec');ylabel('meter/sec^2');
legend('True','Estimate','unCompensated')


figure(2)
subplot(3,1,1)
plot(tFix,yFixErr,'k+-','markersize',5); hold on
plot(t,error.pos,'b-','linewidth',1.5); hold on
% plot(t,xTrue(1,:)-xUnbound(1,:),'g')
plot(t,3*sigma.pos,'r'); hold on
plot(t,3*sigmaMeasureNoise*ones(1,n),'k-.')
plot(t,-3*sigma.pos,'r')
plot(t,-3*sigmaMeasureNoise*ones(1,n),'k-.')
title('Position Error')
xlabel('sec');ylabel('meter');
% legend('Fix Error','Estimate Error','unCompensated Error',...
%     '3\sigma estimate bound')
legend('Fix Error','Estimate Error',...
    '3\sigma estimate bound','3\sigma measurement bound')

% 
subplot(3,1,2)
plot(t,error.vel,'b-','linewidth',1.5); hold on
plot(t,3*sigma.vel,'r'); hold on
plot(t,-3*sigma.vel,'r')
title('Velocity Error')
xlabel('sec');ylabel('meter/sec');
legend('Estimate error','3\sigma estimate bound')


subplot(3,1,3)
plot(t,error.bias,'b-','linewidth',1.5); hold on
plot(t,xTrue(3,:)-xUnbound(3,:),'color',[0 .7 0])
plot(t,3*sigma.bias,'r'); hold on
plot(t,-3*sigma.bias,'r')
legend('Estimate error','Uncompensated error','3\sigma estimate bound')
title('Bias Error')
xlabel('sec');ylabel('meter/sec^2');
%%
% plot of Kalman Gain
figure(3)
subplot(3,1,1)
plot(t,gainStore(1,:),'b','linewidth',1.5); legend('pos gain');
title('Kalman Gains')
grid 
subplot(3,1,2)
plot(t,gainStore(2,:),'b','linewidth',1.5); legend('vel gain'); 
grid 
subplot(3,1,3)
plot(t,gainStore(3,:),'b','linewidth',1.5); legend('bias gain'); 
grid 
