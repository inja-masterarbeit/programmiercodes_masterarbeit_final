function Varianzplot(r, sigma, K, S0, T,  M, M0, N, v0, kappa, theta, rho, N_levels)
format long;
S0 = 100;
K = 100;
T = 1;
M=2;
M0=4;
epsilon=0.01;
V_l=[40.892827019950865, 5.121167149569525, 2.325199583324076, 0.986996681209547, 0.431678536189194, 0.198008343730391, 0.094348272635517, 0.046142300800308, 0.022590165062293, 0.011245767899094, 0.005705397165832, 0.002911580313002]; %asia M=2, M0
%N_l = optimale_Nl(V_l, M, epsilon, T, M0);

r = 0.03; v0 = 0.0225; kappa = 2; theta = 0.0225; sigma = 0.25; rho = -0.5;

N_levels=[2322173 , 424138, 103584,  24122, 5925 , 1465,  369,    91]; 
%gleiche N_l für alle Verfahren für bessere Vergleichbarkeit, direkt sichtbar wer Varianz besser reduziert  

L = length(N_levels) - 1; 
payoffs_diff=zeros(L+1,1);
payoffs_diff_m=zeros(L+1,1);
payoffs_diff_ijk=zeros(L+1,1);

y1=zeros(L+1,1);
y2=zeros(L+1, 1);
y3=zeros(L+1, 1);

payoffs_var=zeros(L+1,1);
payoffs_var_m=zeros(L+1,1);
payoffs_var_ijk=zeros(L+1,1);
payoffs_var_r=zeros(L+1,1);
payoffs_var_m_r=zeros(L+1,1);
 
for l=0:L
    N=N_levels(l+1);
     %N=N_levels(l+1);
     h=T/(M0*M^l);
    
    
     S_fine= zeros(N, M0*M^l+1);     
     S_fine(:,1) = S0;
     v_fine= zeros(N, M0*M^l+1);    
     v_fine(:,1) = v0;
     S_f_m= zeros(N, M0*M^l+1);     
     S_f_m(:,1) = S0;
     v_f_m= zeros(N, M0*M^l+1);    
     v_f_m(:,1) = v0;
     S_f_ijk= zeros(N, M0*M^l+1);     
     S_f_ijk(:,1) = S0;
     v_f_ijk= zeros(N, M0*M^l+1);    
     v_f_ijk(:,1) = v0;
     S_fine_r= zeros(N, M0*M^l+1);     
     S_fine_r(:,1) = S0;
     v_fine_r= zeros(N, M0*M^l+1);    
     v_fine_r(:,1) = v0;
     S_f_m_r= zeros(N, M0*M^l+1);     
     S_f_m_r(:,1) = S0;
     v_f_m_r= zeros(N, M0*M^l+1);    
     v_f_m_r(:,1) = v0;

     Zs=randn(N,M0*M^l);
     Zv=rho*Zs+sqrt(1-rho^2)*randn(N,M0*M^l);  %genriere die beiden Zufallsavriablen mit Cov=rho
     dW1=sqrt(h)*Zs;
     dW2=sqrt(h)*Zv;

 

     for t = 1:M0*M^l
       
        v_fine(:,t+1)=v_fine(:,t)+kappa*(theta-max(0, v_fine(:,t)))*h+sigma*sqrt(max(0, v_fine(:,t))).*dW2(:,t);
        S_fine(:,t+1) = S_fine(:,t) + r * S_fine(:,t) * h + sqrt(max(0, v_fine(:,t))).* S_fine(:,t) .* dW1(:,t);
        
        %Milstein
        v_f_m(:,t+1)=v_f_m(:,t)+kappa*(theta-max(0, v_f_m(:,t)))*h+sigma*sqrt(max(0, v_f_m(:,t))).*dW2(:,t)+1/4*sigma^2*(dW2(:,t).^2-h);
        S_f_m(:,t+1) = S_f_m(:,t) + r * S_f_m(:,t) * h + sqrt(max(0, v_f_m(:,t))).* S_f_m(:,t) .* dW1(:,t)+1/2*S_f_m(:,t).*max(0, v_f_m(:,t)).*(dW1(:,t).^2-h);
        
        %IJK
        v_f_ijk(:,t+1)=(v_f_ijk(:,t)+kappa*theta*h+sigma*sqrt(v_f_ijk(:,t)).*dW2(:,t)+1/4*sigma^2.*(dW2(:,t).^2-h))/(1+kappa*h);
        logS_f_ijk=log(S_f_ijk(:,t))+r*h-1/4*(v_f_ijk(:,t)+v_f_ijk(:,t+1))*h+rho*sqrt(v_f_ijk(:,t)).*dW2(:,t)+1/2*(sqrt(v_f_ijk(:,t))+sqrt(v_f_ijk(:,t+1))).*(dW1(:,t)-rho.*dW2(:,t))+1/4*sigma*rho.*(dW2(:,t).^2-h);
        S_f_ijk(:,t+1) = exp(logS_f_ijk);
         
        %Euler Reflection
        v_fine_r(:,t+1)=v_fine_r(:,t)+kappa*(theta-abs(v_fine_r(:,t)))*h+sigma*sqrt(abs(v_fine_r(:,t))).*dW2(:,t);
        S_fine_r(:,t+1) = S_fine_r(:,t) + r * S_fine_r(:,t) * h + sqrt(abs(v_fine_r(:,t))).* S_fine_r(:,t) .* dW1(:,t);

        %Milstein Reflection
        v_f_m_r(:,t+1)=v_f_m_r(:,t)+kappa*(theta-abs( v_f_m_r(:,t)))*h+sigma*sqrt(abs( v_f_m_r(:,t))).*dW2(:,t)+1/4*sigma^2*(dW2(:,t).^2-h);
        S_f_m_r(:,t+1) = S_f_m_r(:,t) + r * S_f_m_r(:,t) * h + sqrt(abs( v_f_m_r(:,t))).* S_f_m_r(:,t) .* dW1(:,t)+1/2*S_f_m_r(:,t).*abs( v_f_m_r(:,t)).*(dW1(:,t).^2-h);
        
        
 end 

 S_arith_avg = mean(S_fine(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
 payoffs_fine = max(S_arith_avg - K, 0);

 S_arith_avg = mean(S_f_m(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
 payoffs_fine_m = max(S_arith_avg - K, 0);

 S_arith_avg = mean(S_f_ijk(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
 payoffs_fine_ijk = max(S_arith_avg - K, 0);

 S_arith_avg = mean(S_fine_r(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
 payoffs_fine_r = max(S_arith_avg - K, 0);

 S_arith_avg = mean(S_f_m_r(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
 payoffs_fine_m_r = max(S_arith_avg - K, 0);

 y1(l+1)=log(var(payoffs_fine,1))/log(M); %MC Varianzen
 y2(l+1)=log(var(payoffs_fine_m,1))/log(M);
 y3(l+1)=log(var(payoffs_fine_ijk,1))/log(M);

 if l==0
      payoffs_diff(l+1)=mean(payoffs_fine);
      payoffs_diff_m(l+1)=mean(payoffs_fine_m);
      payoffs_diff_ijk(l+1)=mean(payoffs_fine_ijk);
      payoffs_var(l+1)=log(var(payoffs_fine,1))/log(M);
      payoffs_var_m(l+1)=log(var(payoffs_fine_m,1))/log(M);
      payoffs_var_ijk(l+1)=log(var(payoffs_fine_ijk,1))/log(M);
      payoffs_var_r(l+1)=log(var(payoffs_fine_m,1))/log(M);
      payoffs_var_m_r(l+1)=log(var(payoffs_fine_m,1))/log(M);

 else 
      h_coarse=T/(M0*M^(l-1));
      S_coarse = zeros(N,M0* M^(l-1)+1);
      S_coarse(:,1)=S0;
      v_coarse = zeros(N,M0* M^(l-1)+1);
      v_coarse(:,1)=v0;
      S_c_m = zeros(N,M0* M^(l-1)+1);
      S_c_m(:,1)=S0;
      v_c_m = zeros(N,M0* M^(l-1)+1);
      v_c_m(:,1)=v0;
      S_c_ijk = zeros(N,M0* M^(l-1)+1);
      S_c_ijk(:,1)=S0;
      v_c_ijk = zeros(N,M0* M^(l-1)+1);
      v_c_ijk(:,1)=v0;
      S_coarse_r = zeros(N,M0* M^(l-1)+1);
 S_coarse_r(:,1)=S0;
 v_coarse_r = zeros(N,M0* M^(l-1)+1);
 v_coarse_r(:,1)=v0;
  S_c_m_r = zeros(N,M0* M^(l-1)+1);
 S_c_m_r(:,1)=S0;
 v_c_m_r = zeros(N,M0* M^(l-1)+1);
 v_c_m_r(:,1)=v0;
      dW_fine_reshaped = reshape(dW1', M, M0*M^(l-1), N);
      dW1_coarse = squeeze(sum(dW_fine_reshaped, 1))' ;
      dW_fine_reshaped = reshape(dW2', M, M0*M^(l-1), N);
      dW2_coarse = squeeze(sum(dW_fine_reshaped, 1))' ;


     for t=1:M0*M^(l-1)
            v_coarse(:,t+1)=v_coarse(:,t)+kappa*(theta-max(0, v_coarse(:,t)))*h_coarse+sigma*sqrt(max(0, v_coarse(:,t))).*dW2_coarse(:,t);
            S_coarse(:,t+1) = S_coarse(:,t) + r * S_coarse(:,t) * h_coarse + sqrt(max(0, v_coarse(:,t))).* S_coarse(:,t) .* dW1_coarse(:,t);
      
            %Milstein
            v_c_m(:,t+1)=v_c_m(:,t)+kappa*(theta-max(0, v_c_m(:,t)))*h_coarse+sigma*sqrt(max(0, v_c_m(:,t))).*dW2_coarse(:,t)+1/4*sigma^2*(dW2_coarse(:,t).^2-h_coarse);
            S_c_m(:,t+1) = S_c_m(:,t) + r * S_c_m(:,t) * h_coarse + sqrt(max(0, v_c_m(:,t))).* S_c_m(:,t) .* dW1_coarse(:,t)+1/2*S_c_m(:,t).*max(0, v_c_m(:,t)).*(dW1_coarse(:,t).^2-h_coarse);
            
            %IJK
            v_c_ijk(:,t+1)=(v_c_ijk(:,t)+kappa*theta*h_coarse+sigma*sqrt(v_c_ijk(:,t)).*dW2_coarse(:,t)+1/4*sigma^2.*(dW2_coarse(:,t).^2-h_coarse))/(1+kappa*h_coarse);
            logS_c_ijk=log(S_c_ijk(:,t))+r*h_coarse-1/4*(v_c_ijk(:,t)+v_c_ijk(:,t+1))*h_coarse+rho*sqrt(v_c_ijk(:,t)).*dW2_coarse(:,t)+1/2*(sqrt(v_c_ijk(:,t))+sqrt(v_c_ijk(:,t+1))).*(dW1_coarse(:,t)-rho*dW2_coarse(:,t))+1/4*sigma*rho*(dW2_coarse(:,t).^2-h_coarse);
            S_c_ijk(:,t+1)=exp(logS_c_ijk);

            %Euler Reflection
             v_coarse_r(:,t+1)=v_coarse_r(:,t)+kappa*(theta-abs(v_coarse_r(:,t)))*h_coarse+sigma*sqrt(abs(v_coarse_r(:,t))).*dW2_coarse(:,t);
            S_coarse_r(:,t+1) = S_coarse_r(:,t) + r * S_coarse_r(:,t) * h_coarse + sqrt(abs(v_coarse_r(:,t))).* S_coarse_r(:,t) .* dW1_coarse(:,t);

             %Milstein Reflection
            v_c_m_r(:,t+1)=v_c_m_r(:,t)+kappa*(theta-abs(v_c_m_r(:,t)))*h_coarse+sigma*sqrt(abs( v_c_m_r(:,t))).*dW2_coarse(:,t)+1/4*sigma^2*(dW2_coarse(:,t).^2-h_coarse);
            S_c_m_r(:,t+1) = S_c_m_r(:,t) + r * S_c_m_r(:,t) * h_coarse + sqrt(abs(v_c_m_r(:,t))).* S_c_m_r(:,t) .* dW1_coarse(:,t)+1/2*S_c_m_r(:,t).*abs(v_c_m_r(:,t)).*(dW1_coarse(:,t).^2-h_coarse);
            
      
     end

     S_arith_avg = mean(S_coarse(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_coarse = max(S_arith_avg - K, 0);
    
     payoffs_diff(l+1)=mean(payoffs_fine-payoffs_coarse);
     payoffs_var(l+1)=log(var(payoffs_fine-payoffs_coarse,1))/log(M);
    
    
     S_arith_avg = mean(S_c_m(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_coarse_m = max(S_arith_avg - K, 0);

     payoffs_diff_m(l+1)=mean(payoffs_fine_m-payoffs_coarse_m);
     payoffs_var_m(l+1)=log(var(payoffs_fine_m-payoffs_coarse_m,1))/log(M);
     
    
     S_arith_avg = mean(S_c_ijk(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_coarse_ijk = max(S_arith_avg - K, 0);
    
     payoffs_diff_ijk(l+1)=mean(payoffs_fine_ijk-payoffs_coarse_ijk);
     payoffs_var_ijk(l+1)=log(var(payoffs_fine_ijk-payoffs_coarse_ijk,1))/log(M);

     S_arith_avg = mean(S_coarse_r(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_coarse_r = max(S_arith_avg - K, 0);
     payoffs_var_r(l+1)=log(var(payoffs_fine_r-payoffs_coarse_r,1))/log(M);

       S_arith_avg = mean(S_c_m_r(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_coarse_m_r = max(S_arith_avg - K, 0);
     payoffs_var_m_r(l+1)=log(var(payoffs_fine_m_r-payoffs_coarse_m_r,1))/log(M);
     
    
     
 end
 



end
x=0:1:L; 
% Konvergenzrate
P = polyfit(x(2:end), (payoffs_var(2:end)'), 1);
rate = -P(1);              % geschätzte Konvergenzrate (β)
C = exp(P(2));             % geschätzte Konstante C
fprintf('Geschätzte starke Konvergenzordnung Euler: β ≈ %.2f\n', rate);
fprintf('Konstante: C ≈ %.2f\n', C);

P = polyfit(x(2:end), (payoffs_var_m(2:end)'), 1);
rate = -P(1);              % geschätzte Konvergenzrate (β)
C = exp(P(2));             % geschätzte Konstante C
fprintf('Geschätzte starke Konvergenzordnung Milstein: β  ≈ %.2f\n', rate);
fprintf('Konstante: C ≈ %.2f\n', C);

P = polyfit((x(2:end)), (payoffs_var_ijk(2:end)'), 1);
rate = -P(1);              % geschätzte Konvergenzrate (β)
C = exp(P(2));             % geschätzte Konstante C
fprintf('Geschätzte starke Konvergenzordnung Ijk: β  ≈ %.2f\n', rate);
fprintf('Konstante: C ≈ %.2f\n', C);
    
 
figure;
hold on;
colors = lines(5);  % 5 Verfahren × 2 Varianten

% P_l - P_{l-1}
plot(x(2:end), payoffs_var(2:end), '-o', ...
    'Color', colors(1,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'none');
hold on;

plot(x(2:end), payoffs_var_m(2:end), '-o', ...
    'Color', colors(2,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', 'none');

plot(x(2:end), payoffs_var_ijk(2:end), '-o', ...
    'Color', colors(3,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'none');
plot(x(2:end), payoffs_var_r(2:end), '-o', ...
    'Color', colors(4,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(4,:), ...
    'MarkerEdgeColor', 'none');
plot(x(2:end), payoffs_var_m_r(2:end), '-o', ...
    'Color', colors(5,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(5,:), ...
    'MarkerEdgeColor', 'none');

% P_l direkt
plot(x, y1, '-x', 'LineWidth', 2);
plot(x, y2, '-x', 'LineWidth', 2);
plot(x, y3, '-x', 'LineWidth', 2);


hold off;
xlabel('Level l');
ylabel('log_M Varianz');
legend( 'P_l-P_{l-1}(Euler)', 'P_l-P_{l-1}(Milstein)', 'P_l-P_{l-1}(IJK)','P_l-P_{l-1}(Euler, Reflection)','P_l-P_{l-1}(Milstein, Reflection)', ' P_l (Euler)', ' P_l (Milstein)', ' P_l (IJK)');
title('Levelweise Varianzabschätzung (Bestimmung von β)');

end
