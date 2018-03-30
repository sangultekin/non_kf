%Alpha Divergence Kalman Filter
%Algorithm 3 (also Algorithm 2 by setting alpha=1) from the paper:
%S. Gultekin and J. Paisley. Nonlinear Kalman filtering with divergence minimization, 
%IEEE Transactions on Signal Processing, vol. 65, no. 23, pp. 6319-6331.
%Note: Single moment matching step is used for alpha<1 for better accuracy

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com

function [OUT] = func_AKF(IN,MEAS,alpha)

x_filter = IN.x_filter;
P        = IN.P;
F        = IN.F;
Q_CV     = IN.Q_CV;

n = size(x_filter,2);

%number of particles
ns = 1e4;

for t=2:n    
    R = MEAS(t).R;
    data_measurement = MEAS(t).z_true;
    dim = size(x_filter,1);
    
    %predict step
    x_predic = F*x_filter(:,t-1);
    P_predic = F*P*F'+Q_CV;
    
    %use prior as proposal
    x_proposal = x_predic;
    P_proposal = P_predic;
    
    %sample particles
    PART = mvnrnd(x_proposal,P_proposal,ns)';
    PART_sub = PART([1 3],:);
    
    %predict measurements
    tmp_num = MEAS(t).num;
    h = zeros(tmp_num,ns);
    for i=1:tmp_num
        tmp_loc = MEAS(t).loc;
        tmp_dif = PART_sub - repmat(tmp_loc(:,i),1,ns);
        h(i,:) = sqrt( sum(tmp_dif.^2,1) );
    end
    
    %calculate importance weights
    tmp = repmat(data_measurement,1,ns) - h;
    weig = alpha * ( -.5 * sum(tmp.*(inv(R)*tmp) , 1) );
    weig = weig - max(weig);
    weig = exp(weig);
    weig = weig/sum(weig);

    %collapse particles
    mu = PART*weig';
    PART = PART - repmat(mu,1,ns);
    Sig = (PART.*repmat(weig,dim,1))*PART';

    %state estimation
    x_filter(:,t) = mu;
    P = Sig;
end

OUT=x_filter;
end