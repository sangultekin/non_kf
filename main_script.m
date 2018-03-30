%Run all filters and calculate performance.
%Performance comparison of the algorithms presented in the following paper:
%S. Gultekin and J. Paisley. Nonlinear Kalman filtering with divergence minimization, 
%IEEE Transactions on Signal Processing, vol. 65, no. 23, pp. 6319-6331.

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



close all;clear all;clc; 
home=pwd;

%state-space parameters of the data
sample_data_constants;
%overwrite definitions here if you want to simulate parameter uncertainty
%for example the following line feeds a different Q_CV matrix (see paper
%for details)
Q_CV = 1e-1*eye(4); %simply comment out this line if you want to match parameters

%root mean square error over Monte Carlo runs
RMSE_EKF=zeros(1,MC_number);
RMSE_UKF=zeros(1,MC_number);    
RMSE_SIR=zeros(1,MC_number);
RMSE_NKF=zeros(1,MC_number);
RMSE_SKF=zeros(1,MC_number);
RMSE_MKF=zeros(1,MC_number);
RMSE_AKF=zeros(1,MC_number);

%run all filters
for m=1:MC_number    
    cd sample_data;
    load(['data_' num2str(m) '.mat']);
    cd ..;
    
    %prior covariance
    P  = diag([100 .1 100 .1]);
    
    %initialize state estimates
    x_filter_EKF = zeros(4,n);
    x_filter_UKF = zeros(4,n);
    x_filter_SIR = zeros(4,n);
    x_filter_NKF = zeros(4,n);
    x_filter_SKF = zeros(4,n);
    x_filter_MKF = zeros(4,n);
    x_filter_AKF = zeros(4,n);
     
    %fixed starting point to reduce randomness
    t = 1;
    x_filter_EKF(:,t) = target_state - sqrt(diag(P));
    x_filter_UKF(:,t) = x_filter_EKF(:,t);
    x_filter_SIR(:,t) = x_filter_EKF(:,t);
    x_filter_NKF(:,t) = x_filter_EKF(:,t);
    x_filter_SKF(:,t) = x_filter_EKF(:,t);
    x_filter_MKF(:,t) = x_filter_EKF(:,t);
    x_filter_AKF(:,t) = x_filter_EKF(:,t);
    
    %Extended Kalman Filter
    IN_EKF.x_filter = x_filter_EKF;
    IN_EKF.P        = P;
    IN_EKF.F        = F;
    IN_EKF.Q_CV     = Q_CV;
    
    x_filter_EKF = func_EKF(IN_EKF,MEAS);
     
    %Unscented Kalman Filter
    IN_UKF.x_filter = x_filter_UKF;
    IN_UKF.P        = P;
    IN_UKF.F        = F;
    IN_UKF.Q_CV     = Q_CV;
    
    x_filter_UKF = func_UKF(IN_UKF,MEAS);
    
    %Particle Filter: Sequential Importance Resampling version
    IN_SIR.x_filter = x_filter_SIR;
    IN_SIR.P        = P;
    IN_SIR.F        = F;
    IN_SIR.Q_CV     = Q_CV;
    
    x_filter_SIR = func_SIR(IN_SIR,MEAS);
    
    %Ensemble Kalman Filter
    IN_NKF.x_filter = x_filter_NKF;
    IN_NKF.P        = P;
    IN_NKF.F        = F;
    IN_NKF.Q_CV     = Q_CV;
    
    x_filter_NKF = func_NKF(IN_NKF,MEAS);
    
    %Stochastic Search Kalman Filter (Forward KL divergence objective)
    IN_SKF.x_filter = x_filter_SKF;
    IN_SKF.P        = P;
    IN_SKF.F        = F;
    IN_SKF.Q_CV     = Q_CV;
    
    x_filter_SKF = func_SKF(IN_SKF,MEAS);
    
    %Moment Matching Kalman Filter (Reverse KL divergence objective) 
    IN_MKF.x_filter = x_filter_MKF;
    IN_MKF.P        = P;
    IN_MKF.F        = F;
    IN_MKF.Q_CV     = Q_CV;
    
    x_filter_MKF = func_AKF(IN_MKF,MEAS,1); %alpha=1 is MKF
    
    %Alpha Divergence Kalman Filter (Alpha divergence objective)
    IN_AKF.x_filter = x_filter_AKF;
    IN_AKF.P        = P;
    IN_AKF.F        = F;
    IN_AKF.Q_CV     = Q_CV;
    
    x_filter_AKF = func_AKF(IN_AKF,MEAS,.5); %alpha=.5, you can try different values here
    
    %compute errors
    RMSE_EKF(m) = sqrt(  mean( sum( (x_filter_EKF([1 3],:)-data_actual).^2 , 1) )  );
    RMSE_UKF(m) = sqrt(  mean( sum( (x_filter_UKF([1 3],:)-data_actual).^2 , 1) )  );
    RMSE_SIR(m) = sqrt(  mean( sum( (x_filter_SIR([1 3],:)-data_actual).^2 , 1) )  );
    RMSE_NKF(m) = sqrt(  mean( sum( (x_filter_NKF([1 3],:)-data_actual).^2 , 1) )  );
    RMSE_SKF(m) = sqrt(  mean( sum( (x_filter_SKF([1 3],:)-data_actual).^2 , 1) )  ); 
    RMSE_MKF(m) = sqrt(  mean( sum( (x_filter_MKF([1 3],:)-data_actual).^2 , 1) )  ); 
    RMSE_AKF(m) = sqrt(  mean( sum( (x_filter_AKF([1 3],:)-data_actual).^2 , 1) )  );
    
    %print progress
    fprintf('Finished Monte Carlo sample: %i\n', m);
end

%print overall performance
fprintf('\n');
fprintf('Root Mean Square Error Results\n')
fprintf('EKF RMSE: %2.4f\n', mean(RMSE_EKF));
fprintf('UKF RMSE: %2.4f\n', mean(RMSE_UKF));
fprintf('SIR RMSE: %2.4f\n', mean(RMSE_SIR));
fprintf('NKF RMSE: %2.4f\n', mean(RMSE_NKF));
fprintf('SKF RMSE: %2.4f\n', mean(RMSE_SKF));
fprintf('MKF RMSE: %2.4f\n', mean(RMSE_MKF));
fprintf('AKF RMSE: %2.4f\n', mean(RMSE_AKF));