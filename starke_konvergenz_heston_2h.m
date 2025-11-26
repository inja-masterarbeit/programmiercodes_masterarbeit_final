function starke_konvergenz_heston_2h(S0, T, M_list, N, r, v0, kappa, theta, sigma, rho)

% Parameter    
S0 = 100;
T = 1;
r = 0.03; v0 = 0.0225; kappa = 2; theta = 0.0225; sigma = 0.25; rho = -0.5; 
N = 100000;    % Anzahl Pfade
M_list = [8, 16, 32, 64, 128, 256];

errors = zeros(length(M_list), 1);
errors_m = zeros(length(M_list), 1);
errors_ijk = zeros(length(M_list), 1);
for i = 1:length(M_list)
        M = M_list(i);
        h = T / M;


         % Erzeuge gekoppelte Brown'sche Inkremente
        
         Zv=randn(N,M);
         Zs=rho*Zv+sqrt(1-rho^2)*randn(N,M);  %genriere die beiden Zufallsavriablen mit Cov=rho
         dW1=sqrt(h)*Zs;
         dW2=sqrt(h)*Zv;
   

         %Feine Simulation mit h
         S_fine= zeros(N, M+1);     
         S_fine(:,1) = S0;
         v_fine= zeros(N, M+1);    
         v_fine(:,1) = v0;
         S_f_m= zeros(N, M+1);     
         S_f_m(:,1) = S0;
         v_f_m= zeros(N, M+1);    
         v_f_m(:,1) = v0;
         S_f_ijk= zeros(N, M+1);     
         S_f_ijk(:,1) = S0;
         v_f_ijk= zeros(N, M+1);    
         v_f_ijk(:,1) = v0;

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

        % Gröbere Simulation mit 2h
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
                %Euler
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

        % Fehler (Pfadweise Differenz)
        errors(i) = sqrt(mean(abs(S_fine(:,M+1) - S_coarse(:,M/2+1)))^2); %Definition Giles Präsi S. 12 (nur Pfad Ende)
        errors_m(i) = sqrt(mean(abs(S_f_m(:,M+1) - S_c_m(:,M/2+1)))^2); 
        errors_ijk(i) = sqrt(mean(abs(S_f_ijk(:,M+1) - S_c_ijk(:,M/2+1)))^2);
end
% Ausgabe
disp(table(M_list', errors_ijk, 'VariableNames', {'M', 'RMS_Error'}));

% Konvergenzrate
P = polyfit(log(M_list), log(errors'), 1);
rate = -P(1);              % geschätzte Konvergenzrate (β)
C = exp(P(2));             % geschätzte Konstante C
fprintf('Geschätzte starke Konvergenzordnung Euler: β ≈ %.2f\n', rate);
fprintf('Konstante: C ≈ %.2f\n', C);

P = polyfit(log(M_list), log(errors_m'), 1);
rate = -P(1);              % geschätzte Konvergenzrate (β)
C = exp(P(2));             % geschätzte Konstante C
fprintf('Geschätzte starke Konvergenzordnung Milstein: β  ≈ %.2f\n', rate);
fprintf('Konstante: C ≈ %.2f\n', C);

P = polyfit(log(M_list), log(errors_ijk'), 1);
rate = -P(1);              % geschätzte Konvergenzrate (β)
C = exp(P(2));             % geschätzte Konstante C
fprintf('Geschätzte starke Konvergenzordnung Ijk: β  ≈ %.2f\n', rate);
fprintf('Konstante: C ≈ %.2f\n', C);
    
    % Plot
    figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);  % größere Figur
colors = lines(3);  % Farbschema für 3 Verfahren

% Euler
loglog(M_list, errors, 'o-', ...
    'Color', colors(1,:), ...
    'LineWidth', 2, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'none');
hold on;

% Milstein
loglog(M_list, errors_m, 'o-', ...
    'Color', colors(2,:), ...
    'LineWidth', 2, ...
    'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', 'none');
hold on;

% IJK
loglog(M_list, errors_ijk, 'o-', ...
    'Color', colors(3,:), ...
    'LineWidth', 2, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'none');
hold on;
   C = errors(1) * M_list(1)^0.5;   % Konstante so gewählt, dass Linie durch den ersten Punkt geht
   y = 1 *C ./ sqrt(M_list);           % entspricht Steigung -0.5 im loglog-Plot
   loglog(M_list, y, '--o');

    xlabel('Anzahl der Zeitschritte'); ylabel('$\left( E\left[ \left( \hat{S}_{h}(T) - \hat{S}_{2h}(T) \right)^2 \right] \right)^{1/2}$', 'Interpreter', 'latex','FontWeight','bold','FontSize',12);

    title('Starke Konvergenz');
    legend('Euler–Maruyama', 'Milstein', 'IJK', 'Referenzlinie mit Steigung 0.5');
    grid on;

   
end
