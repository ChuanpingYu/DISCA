
clc;
clear;

Num = 50;
R = 1000;
p = 6;
q = 5;
alpha = 0.05;

e = 10^-3;
ninv = norminv(1-alpha/2)^2;

P_U = [0.4082;0.4082;0.4082;0.4082;0.4082;0.4082];
P_V = [0.4472;0.4472;0.4472;0.4472;0.4472];



% mu1 = zeros(p,1);
% sigma1 = [1,0,0;0,1,0.6;0,0.6,1];
% mu2 = zeros(q,1);
% sigma2 = eye(q);
% 
% X_00 = mvnrnd(mu1,sigma1,Num);
% Y_00 = mvnrnd(mu2,sigma2,Num);
% SUM_X = X_00(:,2)+X_00(:,3);
% Y_00(:,1) = 0.01*Y_00(:,1)+SUM_X;

% load /Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/X_1.txt;
% load /Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/Y_1.txt;
% 
% dlmwrite('/Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/one/dimU.txt',[]);
% dlmwrite('/Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/one/dimV.txt',[]);
% dlmwrite('/Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/one/U.txt',[]);
% dlmwrite('/Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/one/V.txt',[]);
% dlmwrite('/Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/one/errorU.txt',[]);
% dlmwrite('/Users/C.Yu/Documents/GeorgiaTech/Research/CCA/Simulations/comparison/one/errorV.txt',[]);

load X50.txt;
load Y50.txt;

dlmwrite('dimU_50_DISCA.txt',[]);
dlmwrite('dimV_50_DISCA.txt',[]);
dlmwrite('U_50_DISCA.txt',[]);
dlmwrite('V_50_DISCA.txt',[]);
dlmwrite('errorU_50_DISCA.txt',[]);
dlmwrite('errorV_50_DISCA.txt',[]);




for r = 1:R
    X_0 = X50((r-1)*Num+(1:Num),:);
    Y_0 = Y50((r-1)*Num+(1:Num),:);

    iter = 0;

    X_diff = zeros(Num*(Num-1)/2,p);  %every row represents one (X_i-X_j) value
    for i = 1:Num
        for j = (i+1):Num
            iter = iter+1;
            X_diff(iter,:) = X_0(i,:)-X_0(j,:);
        end
    end

    U = [];  %each row represents a direction.
    U_n = eye(p);  %null space of U.
    TIME = [];
    DC_X = 0;
    THRES_X = [];

    X_diff_proj = X_diff;
    p_1 = p;

    ind = 0;

    g_0_Y = pdist2(Y_0,Y_0);

    g_Y = zeros(Num*(Num-1)/2,1);  %all g_ij values
    sum3_Y = sum(sum(g_0_Y));
    for i=1:(Num-1)
        sum1 = sum(g_0_Y(i,:));  %sum_k:norm(Y_i-Y_k)
        for j=(i+1):Num
            sum2 = sum(g_0_Y(:,j));  %sum_k:norm(Y_j-Y_k)
            ind = ind +1;
            g_Y(ind) = g_0_Y(i,j)-sum1/Num-sum2/Num+sum3_Y/(Num^2);  %g_ij
        end
    end

    for DIM=1:(p-1)
        X_diff_proj = X_diff*U_n;
        t1 = clock;
        u = cal_min(g_Y,X_diff_proj,Num,p_1,e);
        t2 = clock;
        t = etime(t2,t1);
        u_ult = U_n*u/norm(U_n*u);
        dc_0 = distcov(X_0*u_ult,Y_0);
        DC_X = [DC_X;dc_0];
        p_1 = p_1-1;
        TIME = [TIME;t];
        threshold = ninv*sum3_Y*sum(sum(pdist2(X_0*u_ult,X_0*u_ult)))*(1/Num^5);
        THRES_X = [THRES_X;threshold];
        if dc_0>threshold
            DIM = DIM-1;
            break;
        end
        U = [U;u_ult'];
        U_n=null(U);
        % fprintf('Finished dimension %d of W_1\n',DIM);
    end

    if DIM==p-1
        dc_0 = distcov(X_0*U_n,Y_0);
        threshold = ninv*sum3_Y*sum(sum(pdist2(X_0*U_n,X_0*U_n)))*(1/Num^5);
        DC_X = [DC_X;dc_0];
        THRES_X = [THRES_X;threshold];
        if dc_0<threshold
            DIM = p;
            U_n =[];
        end
    end
    
    dlmwrite('dimU_50_DISCA.txt',p-DIM,'-append');
    dlmwrite('U_50_DISCA.txt',U_n','-append');
    
             

             P1 = U_n*U_n';
             if size(P1,1)==0
             error1=1;
             else
             error1 = norm(P_U*P_U'-P1);
                           end
    dlmwrite('errorU_50_DISCA.txt',error1,'-append');

%     % the optimal dimension of X is 
%     fprintf('The optimal dimension of X is %d.\n',p-DIM);
%     %the basis of the optimal subspace of X is
%     disp('the basis of the optimal subspace of X is (each column is a vector): ');
%     U_n

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iter = 0;

    Y_diff = zeros(Num*(Num-1)/2,q);  %every row represents one (Y_i-Y_j) value
    for i = 1:Num
        for j = (i+1):Num
            iter = iter+1;
            Y_diff(iter,:) = Y_0(i,:)-Y_0(j,:);
        end
    end

    V = [];  %each row represents a direction.
    V_n = eye(q);  %null space of U.
    TIME = [];
    DC_Y = 0;
    THRES_Y = [];

    Y_diff_proj = Y_diff;
    q_1 = q;

    ind = 0;

    if size(U_n,1)==0
        X_p = X_0;
    else
        X_p = X_0*U_n;
    end

    g_0_X = pdist2(X_p,X_p);

    g_X = zeros(Num*(Num-1)/2,1);  %all g_ij values
    sum3_X = sum(sum(g_0_X));
    for i=1:Num
        sum1 = sum(g_0_X(i,:));  %sum_k:norm(Y_i-Y_k)
        for j=(i+1):Num
            sum2 = sum(g_0_X(:,j));  %sum_k:norm(Y_j-Y_k)
            ind = ind +1;
            g_X(ind) = g_0_X(i,j)-sum1/Num-sum2/Num+sum3_X/(Num^2);  %g_ij
        end
    end

    for DIM=1:(q-1)
        Y_diff_proj = Y_diff*V_n;
        t1 = clock;
        v = cal_min(g_X,Y_diff_proj,Num,q_1,e);
        t2 = clock;
        t = etime(t2,t1);
        v_ult = V_n*v/norm(V_n*v);
        dc_0 = distcov(Y_0*v_ult,X_0);
        DC_Y = [DC_Y;dc_0];
        q_1 = q_1-1;
        TIME = [TIME;t];
        threshold = ninv*sum3_X*sum(sum(pdist2(Y_0*v_ult,Y_0*v_ult)))*(1/Num^5);
        THRES_Y = [THRES_Y;threshold];
        if dc_0>threshold
            DIM = DIM-1;
            break;
        end
        V = [V;v_ult'];
        V_n=null(V);
        % fprintf('Finished dimension %d of W_2\n',DIM);
    end

    if DIM==q-1
        dc_0 = distcov(Y_0*V_n,X_0);
        threshold = ninv*sum3_X*sum(sum(pdist2(Y_0*V_n,Y_0*V_n)))*(1/Num^5);
        DC_Y = [DC_Y;dc_0];
        THRES_Y = [THRES_Y;threshold];
        if dc_0<threshold
            DIM = q;
            V_n = [];
        end
    end
    
    dlmwrite('dimV_50_DISCA.txt',q-DIM,'-append');
    dlmwrite('V_50_DISCA.txt',V_n','-append');
    
             P2 = V_n*V_n';
             if size(P2,1)==0
             error2=1;
             else
             error2 = norm(P_V*P_V'-P2);
                           end
    dlmwrite('errorV_50_DISCA.txt',error2,'-append');
    

%     % the optimal dimension of Y is 
%     fprintf('The optimal dimension of Y is %d.\n',q-DIM);
%     %the basis of the optimal subspace of Y is
%     disp('the basis of the optimal subspace of Y is (each column is a vector): ');
%     V_n

end


