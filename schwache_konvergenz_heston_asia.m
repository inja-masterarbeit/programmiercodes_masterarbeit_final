function schwache_konvergenz_heston_asia
    % Parameter
    
S0 = 100;
K = 100;
T=1;
P=5;
M_list = [8, 16, 32, 64, 128, 254];


r = 0.03; v0 = 0.0225; kappa = 2; theta = 0.0225; sigma = 0.25; rho = -0.5;
%sigma = 0.9;  rho = -0.9;  kappa = 0.5; %stärkere Koeff für mehr neg
errors = zeros(length(M_list), P);
        errors_m = zeros(length(M_list), P);
        errors_ijk = zeros(length(M_list), P);
    for p=1:P
        N = 100000;                         % Anzahl Pfade
        
       
        for i = 1:length(M_list)
            M = M_list(i);
            h = T / M;


            % Erzeuge gekoppelte Brown'sche Inkremente
            
            Zv=randn(N,M);
            Zs=rho*Zv+sqrt(1-rho^2)*randn(N,M);  %genriere die beiden Zufallsavriablen mit Cov=rho
            dW1=sqrt(h)*Zs;
            dW2=sqrt(h)*Zv;
   

            %Feine Simulation mit h
            S_fine=zeros(N,M+1);
            S_fine(:,1)=S0;
            v_fine=zeros(N,M+1);
            v_fine(:,1)=v0;
            S_f_m=zeros(N,M+1);
            S_f_m(:,1)=S0;
            v_f_m=zeros(N,M+1);
            v_f_m(:,1)=v0;
            S_f_ijk=zeros(N,M+1);
            S_f_ijk(:,1)=S0;
            v_f_ijk=zeros(N,M+1);
            v_f_ijk(:,1)=v0;
            for t = 1:M
                    v_fine(:,t+1)=v_fine(:,t)+kappa*(theta-max(0, v_fine(:,t)))*h+sigma*sqrt(max(0, v_fine(:,t))).*dW2(:,t);
                    S_fine(:,t+1) = S_fine(:,t) + r * S_fine(:,t) * h + sqrt(max(0, v_fine(:,t))).* S_fine(:,t) .* dW1(:,t);
                    
                    %Milstein
                    v_f_m(:,t+1)=v_f_m(:,t)+kappa*(theta-max(0, v_f_m(:,t)))*h+sigma*sqrt(max(0, v_f_m(:,t))).*dW2(:,t)+1/4*sigma^2*(dW2(:,t).^2-h);
                    S_f_m(:,t+1) = S_f_m(:,t) + r * S_f_m(:,t) * h + sqrt(max(0, v_f_m(:,t))).* S_f_m(:,t) .* dW1(:,t)+1/2*S_f_m(:,t).*max(0, v_f_m(:,t)).*(dW1(:,t).^2-h);
                    
                    %IJK
                    v_f_ijk(:,t+1)=(v_f_ijk(:,t)+kappa*theta*h+sigma*sqrt(v_f_ijk(:,t)).*dW2(:,t)+1/4*sigma^2.*(dW2(:,t).^2-h))/(1+kappa*h);
                    logS_f_ijk=log(S_f_ijk(:,t))+r*h-1/4*(v_f_ijk(:,t)+v_f_ijk(:,t+1))*h+rho*sqrt(v_f_ijk(:,t)).*dW2(:,t)+1/2*(sqrt(v_f_ijk(:,t))+sqrt(v_f_ijk(:,t+1))).*(dW1(:,t)-rho.*dW2(:,t))+1/4*sigma*rho.*(dW2(:,t).^2-h);
                    S_f_ijk(:,t+1) = exp(logS_f_ijk);
            end 
             S_arith_avg = mean(S_fine(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_fine = max(S_arith_avg - K, 0);
    
     S_arith_avg = mean(S_f_m(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_f_m = max(S_arith_avg - K, 0);

     
     
    
     S_arith_avg = mean(S_f_ijk(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_f_ijk = max(S_arith_avg - K, 0);
        
            % Grobe Simulation mit 2h
            h_coarse = 2*h;
            dW1_coarse = dW1(:,1:2:end) + dW1(:,2:2:end);  % Paarweise summiert um aus h -> 2h, gekoppelten Brownschen Inkremente
            dW2_coarse = dW2(:,1:2:end) + dW2(:,2:2:end);
            
            S_coarse=zeros(N,M/2+1);
            S_coarse(:,1)=S0;
            v_coarse=zeros(N,M/2+1);
            v_coarse(:,1)=v0;
            S_c_m=zeros(N,M/2+1);
            S_c_m(:,1)=S0;
            v_c_m=zeros(N,M/2+1);
            v_c_m(:,1)=v0;
            S_c_ijk=zeros(N,M/2+1);
            S_c_ijk(:,1)=S0;
            v_c_ijk=zeros(N,M/2+1);
            v_c_ijk(:,1)=v0;
            for t = 1:(M/2)
                    v_coarse(:,t+1)=v_coarse(:,t)+kappa*(theta-max(0, v_coarse(:,t)))*h_coarse+sigma*sqrt(max(0, v_coarse(:,t))).*dW2_coarse(:,t);
                    S_coarse(:,t+1) = S_coarse(:,t) + r * S_coarse(:,t) * h_coarse + sqrt(max(0, v_coarse(:,t))).* S_coarse(:,t) .* dW1_coarse(:,t);
              
                    %Milstein
                    v_c_m(:,t+1)=v_c_m(:,t)+kappa*(theta-max(0, v_c_m(:,t)))*h_coarse+sigma*sqrt(max(0, v_c_m(:,t))).*dW2_coarse(:,t)+1/4*sigma^2*(dW2_coarse(:,t).^2-h_coarse);
                    S_c_m(:,t+1) = S_c_m(:,t) + r * S_c_m(:,t) * h_coarse + sqrt(max(0, v_c_m(:,t))).* S_c_m(:,t) .* dW1_coarse(:,t)+1/2*S_c_m(:,t).*max(0, v_c_m(:,t)).*(dW1_coarse(:,t).^2-h_coarse);
                    
                    %IJK
                    v_c_ijk(:,t+1)=(v_c_ijk(:,t)+kappa*theta*h_coarse+sigma*sqrt(v_c_ijk(:,t)).*dW2_coarse(:,t)+1/4*sigma^2.*(dW2_coarse(:,t).^2-h_coarse))/(1+kappa*h_coarse);
                    logS_c_ijk=log(S_c_ijk(:,t))+r*h_coarse-1/4*(v_c_ijk(:,t)+v_c_ijk(:,t+1))*h_coarse+rho*sqrt(v_c_ijk(:,t)).*dW2_coarse(:,t)+1/2*(sqrt(v_c_ijk(:,t))+sqrt(v_c_ijk(:,t+1))).*(dW1_coarse(:,t)-rho*dW2_coarse(:,t))+1/4*sigma*rho*(dW2_coarse(:,t).^2-h_coarse);
                    S_c_ijk(:,t+1)=exp(logS_c_ijk);
            end
           
             S_arith_avg = mean(S_coarse(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_coarse = max(S_arith_avg - K, 0);
    
     S_arith_avg = mean(S_c_m(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_c_m = max(S_arith_avg - K, 0);

     
     
    
     S_arith_avg = mean(S_c_ijk(:, 2:end), 2);  % Arithmetisches Mittel über die Zeit
     payoffs_c_ijk = max(S_arith_avg - K, 0);
         

            errors(i, p) = abs(mean(payoffs_fine)-mean(payoffs_coarse));
            errors_m(i,p) = abs(mean(payoffs_f_m) -mean( payoffs_c_m));
            errors_ijk(i,p) = abs(mean(payoffs_f_ijk )-mean( payoffs_c_ijk));


   
        end

    end
    errors_ges=sum(errors,2)*1/P; %Durschnitt über die P Durchläufe, um schwache Konvergenz zu erhalten
    errors_m_ges=sum(errors_m,2)*1/P;
    errors_ijk_ges=sum(errors_ijk,2)*1/P;
    % Ausgabe
    disp(table(M_list', errors_ges, 'VariableNames', {'M', 'Error Euler'}));
    disp(table(M_list', errors_m_ges, 'VariableNames', {'M', 'Error Milstein'}));
    disp(table(M_list', errors_ijk_ges, 'VariableNames', {'M', 'Error IJK'}));
 

   
    figure;
    loglog(M_list, errors_ges, 'o-', 'LineWidth', 2);
    hold on;
    
    loglog(M_list, errors_m_ges, 'o-', 'LineWidth', 2);
    hold on;
    
    loglog(M_list, errors_ijk_ges, 'o-', 'LineWidth', 2);
    hold on;

    C = errors_ges(1)*M_list(1);   
    y = 1.5*C ./ M_list;               
    loglog(M_list, y, '--o');

   
    xlabel('Anzahl Zeitschritte M'); ylabel('$\mid E\left[ \left( f(\hat{S}_{h}) -f( \hat{S}_{2h}) \right)\right] \mid$', 'Interpreter', 'latex','FontWeight','bold','FontSize',12);

    title('Schwache Konvergenz ');
    legend('Euler', 'Milstein', 'IJK', 'Referenzlinie mit Steigung 1');
    grid on;

end