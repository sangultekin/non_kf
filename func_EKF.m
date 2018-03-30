%Extended Kalman Filter

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function OUT = func_EKF(IN,MEAS)

x_filter = IN.x_filter;
P        = IN.P;
F        = IN.F;
Q_CV     = IN.Q_CV;

n = size(x_filter,2);

%t=1 is initial position
for t=2:n
    R = MEAS(t).R;
    data_measurement = MEAS(t).z_true;
    
    x_predic = F*x_filter(:,t-1);
    P_predic = F*P*F'+Q_CV;
    
    %linearize with auxiliary jacobian script
    [H, z_predic] = aux_jac(x_predic,MEAS(t).loc);
    
    S=H*P_predic*H'+R;
    K=P_predic*H'*inv(S);
    x_filter(:,t)=x_predic+K*(data_measurement-z_predic);
    P=P_predic-K*S*K';
end

OUT=x_filter;
end