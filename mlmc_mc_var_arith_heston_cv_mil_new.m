function price1 = mlmc_mc_var_arith_heston_cv_mil_new (r, sigma, K, S0, T,h,  M, N, dW)
%arithmetische im Heston, Geometrische im Black-Scholes
S0 = 100;
K = 100;
T = 1;
M=2;
M0=4;
 
epsilon=0.01;
 
r = 0.03; v0 = 0.0225; kappa = 2; theta = 0.0225; sigma = 0.25; rho = -0.5; 
V_l=[40.892827019950865, 5.121167149569525, 2.325199583324076, 0.986996681209547, 0.431678536189194, 0.198008343730391, 0.094348272635517, 0.046142300800308, 0.022590165062293, 0.011245767899094, 0.005705397165832, 0.002911580313002]; %asia M=2, M0
N_l = optimale_Nl(V_l, M, epsilon, T, M0);
L = 10; 

price_cv=0;
var_cv_fine=zeros(L+1,1);
var_mlmc_cv=zeros(L,1);
payoffs_var_m=zeros(L+1,1);
for l=0:L
     N=N_l(l+1);
     S_fine = zeros(N,M0* M^l+1);         % Speicher für Pfad
     v_fine = zeros(N,M0*M^l+1);
     v_fine(:,1)=v0;
     S_fine(:,1) = S0;     % Startwert für jeden Pfad
     S_f_bs=zeros(N,M0*M^l+1);
     S_f_bs(:,1)=S0;
     h=T/(M0*M^l);
    
     n=M0*M^l;

     %% Analytischer Preis (diskret, Kemna & Vorst)
    
    sigma_eq = sqrt(theta + (v0 - theta) / (kappa * T) * (1 - exp(-kappa * T)));
    m  = log(S0) + (r - 0.5*sigma_eq^2) * ((n+1)/(2*n)) * T;
    s2 = sigma_eq^2 * T * ((n+1)*(2*n+1)) / (6*n^2);
    s  = sqrt(s2);
    
    d1 = (m + s2 - log(K)) / s;
    d2 = (m - log(K)) / s;
    
    price_exact = (exp(m + 0.5*s2) * normcdf(d1) - K * normcdf(d2));
    
    Zs=randn(N,M0*M^l);
    Zv=rho*Zs+sqrt(1-rho^2)*randn(N,M0*M^l);  %genriere die beiden Zufallsavriablen mit Cov=rho
    dW1=sqrt(h)*Zs;
    dW2=sqrt(h)*Zv;

     for t = 1:M0*M^l
            v_fine(:,t+1)=v_fine(:,t)+kappa*(theta-max(0, v_fine(:,t)))*h+sigma*sqrt(max(0, v_fine(:,t))).*dW2(:,t)+1/4*sigma^2*(dW2(:,t).^2-h);
            S_fine(:,t+1) = S_fine(:,t) + r * S_fine(:,t) * h + sqrt(max(0, v_fine(:,t))).* S_fine(:,t) .* dW1(:,t)+1/2*S_fine(:,t).*max(0, v_fine(:,t)).*(dW1(:,t).^2-h);
            
            S_f_bs(:,t+1) = S_f_bs(:,t) .* exp((r - 0.5*sigma_eq^2)*h + sigma_eq.*dW1(:,t)); %exakte Brownsche Bewegung
     end
    % Beobachtungspunkte gleichmäßig verteilt
    idx = floor(linspace(1, M0*M^l+1, n+1)); 
    idx = idx(2:end);
    S_obs = S_f_bs(:, idx);
    % Geometrisches Mittel
    S_geom = exp(mean(log(S_obs), 2));
    payoffs_MC = max(S_geom - K, 0);
    
    S_arith_avg = mean(S_fine(:, 2:end), 2);
    payoffs_arith= max(S_arith_avg - K, 0);
    b_l = cov(payoffs_arith, payoffs_MC);
    b_l = b_l(1,2) / var(payoffs_MC);
    if l==0
            price_cv=exp(-r*T)*mean(payoffs_arith-b_l*(payoffs_MC-price_exact));
            
            payoffs_var_m(l+1)=log(var(payoffs_arith)); %Var MLMC
    
            var_total_mlmc_cv=var(payoffs_arith-b_l*(payoffs_MC-price_exact))/N; %Gesamtvarian MLMC+CV
            var_cv_fine(l+1)=log(var(payoffs_arith-b_l*(payoffs_MC-price_exact))); 
            var_mc=var((payoffs_arith))/N; %Gesamtvarianz MLMC
     else
             S_coarse = zeros(N, M0*M^(l-1)+1);
             v_coarse=zeros(N,M0*M^(l-1)+1);
             S_coarse(:,1)=S0;
             v_coarse(:,1)=v0;
             S_c_bs = zeros(N, M0*M^(l-1)+1);
             S_c_bs(:,1)=S0;
             h_coarse=T/(M0*M^(l-1));
    
             n=M0*M^(l-1);
            %% Analytischer Preis (diskret, Kemna & Vorst)
            
            sigma_eq = sqrt(theta + (v0 - theta) / (kappa * T) * (1 - exp(-kappa * T)));
            m  = log(S0) + (r - 0.5*sigma_eq^2) * ((n+1)/(2*n)) * T;
            s2 = sigma_eq^2 * T * ((n+1)*(2*n+1)) / (6*n^2);
            s  = sqrt(s2);
            
            d1 = (m + s2 - log(K)) / s;
            d2 = (m - log(K)) / s;
            
            price_exact_grob =  (exp(m + 0.5*s2) * normcdf(d1) - K * normcdf(d2));
            dW1_fine_reshaped = reshape(dW1', M,M0* M^(l-1), N);
            dW2_fine_reshaped = reshape(dW2', M,M0* M^(l-1), N);
            dW1_coarse = squeeze(sum(dW1_fine_reshaped, 1))'; 
            dW2_coarse = squeeze(sum(dW2_fine_reshaped, 1))'; 
    
             for t=1:M0*M^(l-1)
                   v_coarse(:,t+1)=v_coarse(:,t)+kappa*(theta-max(0, v_coarse(:,t)))*h_coarse+sigma*sqrt(max(0, v_coarse(:,t))).*dW2_coarse(:,t)+1/4*sigma^2*(dW2_coarse(:,t).^2-h_coarse);
                   S_coarse(:,t+1) = S_coarse(:,t) + r * S_coarse(:,t) * h_coarse + sqrt(max(0, v_coarse(:,t))).* S_coarse(:,t) .* dW1_coarse(:,t)+1/2*S_coarse(:,t).*max(0, v_coarse(:,t)).*(dW1_coarse(:,t).^2-h_coarse);
                
                   S_c_bs(:,t+1) = S_c_bs(:,t) .* exp((r - 0.5*sigma_eq^2)*h_coarse + sigma_eq.*dW1_coarse(:,t));
             end
            % Beobachtungspunkte gleichmäßig über [T/n, ..., T]
            idx = floor(linspace(1, M0*M^(l-1)+1, n+1)); 
            idx = idx(2:end);
            S_obs = S_c_bs(:, idx);
            % Geometrisches Mittel
            S_geom = exp(mean(log(S_obs), 2));
            payoffs_c_MC = max(S_geom - K, 0);
    
            S_arith_avg = mean(S_coarse(:, 2:end), 2);
            payoffs_c_arith= max(S_arith_avg - K, 0);
    
            b_l = cov(payoffs_arith-payoffs_c_arith, payoffs_MC- payoffs_c_MC);
            b_l = b_l(1,2) / var(payoffs_MC- payoffs_c_MC);
            b_l_f = cov(payoffs_arith, payoffs_MC);
            b_l_f = b_l_f(1,2) / var(payoffs_MC);
           
            price_cv=price_cv+exp(-r*T)*mean((payoffs_arith-payoffs_c_arith)-b_l*(payoffs_MC- payoffs_c_MC - (price_exact-price_exact_grob)));
    
    
            var_mlmc_cv(l)=log(var((payoffs_arith-payoffs_c_arith)-b_l*(payoffs_MC- payoffs_c_MC-(price_exact-price_exact_grob)))); %Var MLMC + CV
            
            var_cv_fine(l+1)=log(var((payoffs_arith)-b_l_f*(payoffs_MC-(price_exact)))); %Var CV
            payoffs_var_m(l+1)=log(var(payoffs_arith-payoffs_c_arith)); %Var MLMC
             
            var_total_mlmc_cv=var_total_mlmc_cv+var((payoffs_arith-payoffs_c_arith)-b_l*(payoffs_MC- payoffs_c_MC-(price_exact-price_exact_grob)))/N; %Gesamtvarianz MLMC+CV
            var_mc=var_mc+var((payoffs_arith-payoffs_c_arith))/N; %Gesamtvarianz MLMC
       end
end
price1=price_cv;


payoffs_diff=zeros(L+1,1); 
payoffs_diff_ijk=zeros(L+1,1);
y1=zeros(L+1,1);
y2=zeros(L+1, 1);
y3=zeros(L+1, 1);
payoffs_var=zeros(L+1,1);
payoffs_var_ijk=zeros(L+1,1);
 
for l=0:L
    N=N_l(l+1);
    h=T/(M0*M^l);


     S_fine= zeros(N, M0*M^l+1);     
     S_fine(:,1) = S0;
     v_fine= zeros(N, M0*M^l+1);    
     v_fine(:,1) = v0;
     S_f_ijk= zeros(N, M0*M^l+1);     
     S_f_ijk(:,1) = S0;
     v_f_ijk= zeros(N, M0*M^l+1);    
     v_f_ijk(:,1) = v0;

     Zs=randn(N,M0*M^l);
     Zv=rho*Zs+sqrt(1-rho^2)*randn(N,M0*M^l);  %genriere die beiden Zufallsavriablen mit Cov=rho
     dW1=sqrt(h)*Zs;
     dW2=sqrt(h)*Zv;

 

     for t = 1:M0*M^l
      
            v_fine(:,t+1)=v_fine(:,t)+kappa*(theta-max(0, v_fine(:,t)))*h+sigma*sqrt(max(0, v_fine(:,t))).*dW2(:,t);
            S_fine(:,t+1) = S_fine(:,t) + r * S_fine(:,t) * h + sqrt(max(0, v_fine(:,t))).* S_fine(:,t) .* dW1(:,t);
            
            %IJK
            v_f_ijk(:,t+1)=(v_f_ijk(:,t)+kappa*theta*h+sigma*sqrt(v_f_ijk(:,t)).*dW2(:,t)+1/4*sigma^2.*(dW2(:,t).^2-h))/(1+kappa*h);
            logS_f_ijk=log(S_f_ijk(:,t))+r*h-1/4*(v_f_ijk(:,t)+v_f_ijk(:,t+1))*h+rho*sqrt(v_f_ijk(:,t)).*dW2(:,t)+1/2*(sqrt(v_f_ijk(:,t))+sqrt(v_f_ijk(:,t+1))).*(dW1(:,t)-rho.*dW2(:,t))+1/4*sigma*rho.*(dW2(:,t).^2-h);
            S_f_ijk(:,t+1) = exp(logS_f_ijk);
     end 
    
    S_arith_avg = mean(S_fine(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
    payoffs_fine = max(S_arith_avg - K, 0);
    
    S_arith_avg = mean(S_f_ijk(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
    payoffs_fine_ijk = max(S_arith_avg - K, 0);
    
    y1(l+1)=log(var(payoffs_fine));
    y3(l+1)=log(var(payoffs_fine_ijk));
    
     if l==0
            payoffs_diff(l+1)=mean(payoffs_fine);
            payoffs_var(l+1)=log(var(payoffs_fine));
          
     else 
             h_coarse=T/(M0*M^(l-1));
            
             S_coarse = zeros(N,M0* M^(l-1)+1);
             S_coarse(:,1)=S0;
             v_coarse = zeros(N,M0* M^(l-1)+1);
             v_coarse(:,1)=v0;
             S_c_ijk = zeros(N,M0* M^(l-1)+1);
             S_c_ijk(:,1)=S0;
             v_c_ijk = zeros(N,M0* M^(l-1)+1);
             v_c_ijk(:,1)=v0;
            
             dW_fine_reshaped = reshape(dW1', M, M0*M^(l-1), N);
             dW1_coarse = squeeze(sum(dW_fine_reshaped, 1))' ;
             dW_fine_reshaped = reshape(dW2', M, M0*M^(l-1), N);
             dW2_coarse = squeeze(sum(dW_fine_reshaped, 1))' ;
            
            
            
             for t=1:M0*M^(l-1)
                        v_coarse(:,t+1)=v_coarse(:,t)+kappa*(theta-max(0, v_coarse(:,t)))*h_coarse+sigma*sqrt(max(0, v_coarse(:,t))).*dW2_coarse(:,t);
                        S_coarse(:,t+1) = S_coarse(:,t) + r * S_coarse(:,t) * h_coarse + sqrt(max(0, v_coarse(:,t))).* S_coarse(:,t) .* dW1_coarse(:,t);
                  
                        %IJK
                        v_c_ijk(:,t+1)=(v_c_ijk(:,t)+kappa*theta*h_coarse+sigma*sqrt(v_c_ijk(:,t)).*dW2_coarse(:,t)+1/4*sigma^2.*(dW2_coarse(:,t).^2-h_coarse))/(1+kappa*h_coarse);
                        logS_c_ijk=log(S_c_ijk(:,t))+r*h_coarse-1/4*(v_c_ijk(:,t)+v_c_ijk(:,t+1))*h_coarse+rho*sqrt(v_c_ijk(:,t)).*dW2_coarse(:,t)+1/2*(sqrt(v_c_ijk(:,t))+sqrt(v_c_ijk(:,t+1))).*(dW1_coarse(:,t)-rho*dW2_coarse(:,t))+1/4*sigma*rho*(dW2_coarse(:,t).^2-h_coarse);
                        S_c_ijk(:,t+1)=exp(logS_c_ijk);
             end
            
             S_arith_avg = mean(S_coarse(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
             payoffs_coarse = max(S_arith_avg - K, 0);
            
             payoffs_diff(l+1)=mean(payoffs_fine-payoffs_coarse);
             payoffs_var(l+1)=log(var(payoffs_fine-payoffs_coarse));
            
            
            S_arith_avg = mean(S_c_ijk(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
            payoffs_coarse_ijk = max(S_arith_avg - K, 0);
            
            payoffs_diff_ijk(l+1)=mean(payoffs_fine_ijk-payoffs_coarse_ijk);
            payoffs_var_ijk(l+1)=log(var(payoffs_fine_ijk-payoffs_coarse_ijk));
     end
     



end
figure;
x=0:1:L;
disp(size(x(2:end)));

plot(x(2:end), payoffs_var(2:end),'-o',  'LineWidth', 1);
hold on;
plot(x(2:end), payoffs_var_m(2:end), '-o',  'LineWidth', 1);
hold on;
plot(x(2:end), payoffs_var_ijk(2:end), '-o',  'LineWidth', 1) ;
hold on;
plot(x(2:end), var_mlmc_cv, '-o', 'Color', 'black',  'LineWidth', 1);
hold on;
plot(x, y1, '-x',  'LineWidth', 1);
hold on;

plot(x, var_cv_fine, '-x',  'LineWidth', 1);
xlabel('Level l');
ylabel('log_M Varianz');
legend( 'P^l-P^{l-1}(Euler)',  'P^l-P^{l-1}(Milstein)', 'P^l-P^{l-1}(IJK)', 'P^l-P^{l-1} (Euler, Cv)', 'P^l (Euler)', 'P^l (Euler, CV)');
title('Levelweise Varianzabschätzung mit CV (Euler)');



