%Unscented Kalman Filter

%Author: San Gultekin, Columbia University, san.gultekin@gmail.com



function OUT = func_UKF(IN,MEAS)

x_filter = IN.x_filter;
P        = IN.P;
F        = IN.F;
Q_CV     = IN.Q_CV;

n = size(x_filter,2);

%UKF parameters
npt=4; 
lamda=2-npt;
Wm = repmat(1/(2*(npt+lamda)), 1,2*npt+1);
Wc = repmat(1/(2*(npt+lamda)), 1,2*npt+1);
Wm(1)=lamda/(npt+lamda);
Wc(1)=lamda/(npt+lamda);

for t=2:n
    R = MEAS(t).R;
    data_measurement = MEAS(t).z_true;
    dim = size(x_filter,1);
    dim2 = MEAS(t).num;
    
    %generate sigma points
    P = .5*(P+P'); cho1 = chol( P )';
    xgamaP1 = repmat(x_filter(:,t-1),1,npt) + sqrt(npt+lamda)*cho1;
    xgamaP2 = repmat(x_filter(:,t-1),1,npt) - sqrt(npt+lamda)*cho1;
    Xsigma=[x_filter(:,t-1) xgamaP1 xgamaP2];

    %predict step
    Xsigma_predic=F*Xsigma;
    x_predic = sum( bsxfun(@times,Wm,Xsigma_predic) , 2);
    P_predic=zeros(size(P));
    for j=1:2*npt+1
        P_predic=P_predic+Wc(j)*( (Xsigma_predic(:,j)-x_predic)*(Xsigma_predic(:,j)-x_predic)' );
    end
    P_predic=P_predic+Q_CV;   
    
    %predict measurements from sigma points
    Zsigma_predic = zeros(dim2,2*npt+1);
    for j=1:2*npt+1
        tmp_vec = Xsigma_predic([1 3],j);
        tmp_dif = repmat(tmp_vec,1,dim2) - MEAS(t).loc;
        Zsigma_predic(:,j)=sqrt(sum(tmp_dif.^2,1))';
    end
    z_predic = sum( bsxfun(@times,Wm,Zsigma_predic) , 2);
 
    %update step
    Pzz=zeros(dim2,dim2);
    for j=1:2*npt+1
        Pzz=Pzz+Wc(j)*( (Zsigma_predic(:,j)-z_predic)*(Zsigma_predic(:,j)-z_predic)' );
    end
    Pzz=Pzz+R;
    
    Pxz=zeros(dim,dim2);
    for j=1:2*npt+1
        Pxz=Pxz+Wc(j)*(Xsigma_predic(:,j)-x_predic)*(Zsigma_predic(:,j)-z_predic)';
    end
    
    K=Pxz*inv(Pzz);
    x_filter(:,t)=x_predic+K*(data_measurement-z_predic);
    P=P_predic-K*Pzz*K';
end

OUT=x_filter;
end