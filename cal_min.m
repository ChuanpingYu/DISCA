function u_pre = cal_min(g,X_diff,Num,p,e)

    xi = 200;
    psi = 10; 
    omega = 10;
    epsilon_abs = 10^-3;
    epsilon_rel = 10^-3;    
        
    ind_pos = find(g>0);
    ind_neg = find(g<0);
    M_plus = repmat(g(ind_pos),[1,p]).*X_diff(ind_pos,:);
    ind_zero_M_plus = find(any(M_plus,2)==0);
    M_plus(ind_zero_M_plus,:) = [];
    M_minus = repmat(-g(ind_neg),[1,p]).*X_diff(ind_neg,:);
    ind_zero_M_minus = find(any(M_minus,2)==0);
    M_minus(ind_zero_M_minus,:) = [];

    %u_pre = (X_diff'*g)/norm(X_diff'*g);  %initial u 
    u_pre = zeros(p,1);
    u_pre(1) = 1;
    
    if max(isnan(u_pre))
        u_pre = zeros(p,1);
        for w=1:p
            s = num2str(X_diff(:,w)'*g,'%1.4e');
            id = find(s=='e');
            s = s(1:id-1);
            u_pre(w) = str2double(s);
        end
        u_pre = u_pre/norm(u_pre);
    end
    
    u = zeros(p,1);
    index = 0;

    while min(norm(u-u_pre),norm(u+u_pre))/max(norm(u),1)>e
        index = index + 1;
        %display(index);
        u = u_pre;
        sig = sign(M_minus*u);
        partial = sig;
        sig_zero = sig(sig==0);
        if sum(sig==0)>0
            sig_zero = 2*rand(size(sig_zero,1),1)-1;
            partial(sig==0) = sig_zero;
        end
        if u==zeros(p,1)
            y = M_minus'*partial;
        else
            y = ((xi-psi)/norm(u))*u+M_minus'*partial;
        end
        rho = 1000;
        z = M_plus*u;
        v = 10*ones(size(M_plus,1),1);
        for l=1:10^5
            u_l = (xi*eye(p)+rho*(M_plus'*M_plus))\(y+M_plus'*(rho*z-v));
            x = (1/rho)*v+M_plus*u_l;
            z_l = sign(x).*max(abs(x)-1/rho,0);
            v = v+rho*(M_plus*u_l-z_l);
            r = M_plus*u_l-z_l;
            s = rho*M_plus'*(z_l-z);
            if norm(r)<=sqrt(size(M_plus,1))*epsilon_abs+epsilon_rel*max(norm(M_plus*u_l),norm(z_l)) &&...
                    norm(s)<=sqrt(p)*epsilon_abs+epsilon_rel*norm(M_plus'*v)
                break
            end
            z = z_l;
        end
        u_pre = u_l;
        psi = psi + xi*(norm(u_pre)-1);
        xi = omega*xi;
    end