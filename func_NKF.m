%Ensemble Kalman Filter

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function [OUT] = func_NKF(IN,MEAS)

x_filter = IN.x_filter;
P        = IN.P;
F        = IN.F;
Q_CV     = IN.Q_CV;

n = size(x_filter,2);

ns = 1e4;
PART = mvnrnd(x_filter(:,1),P,ns)';

for t=2:n
    R = MEAS(t).R;
    data_measurement = MEAS(t).z_true;
    dim = size(x_filter,1);
    
    %propagate particles
    PART = F*PART + mvnrnd(zeros(dim,1),Q_CV,ns)';
    PART_sub = PART([1 3],:);
    
    %predict measurements
    tmp_num = MEAS(t).num;
    h = zeros(tmp_num,ns);
    for i=1:tmp_num
        tmp_loc = MEAS(t).loc;
        tmp_dif = PART_sub - repmat(tmp_loc(:,i),1,ns);
        h(i,:) = sqrt( sum(tmp_dif.^2,1) );
    end
    
    %data assimilation
    PART_Y = h + sqrt(R)*randn(size(h));
    
    REP_Y = repmat(data_measurement,1,ns);
    
    EXf = PART - repmat( sum(PART,2) , 1 , size(PART,2) )/ns;
    EYf = PART_Y - repmat( sum(PART_Y,2) , 1 , size(PART_Y,2) )/ns;
    
    PXX = (EXf * EXf') / (ns-1);
    PYY = (EYf * EYf') / (ns-1);
    PXY = (EXf * EYf') / (ns-1);
    
    K = PXY * inv(PYY);  
	
    PART = PART + K*(REP_Y - PART_Y);
    
    %state estimation
    x_filter(:,t) = sum(PART , 2) / ns;
end

OUT=x_filter;
end