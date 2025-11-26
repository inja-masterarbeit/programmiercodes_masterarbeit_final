function N_opt = Nl_vs_l(V_l, V_l_m, V_l_i, M, M0, epsilon, T, K, r, v0, kappa, theta, sigma, rho, S0)

K = 100; T = 1.0; 
r = 0.03; v0 = 0.0225; kappa = 2; theta = 0.0225;sigma = 0.25;rho = -0.5; 

M=2; M0=4; S0 = 100;
epsilon=0.01;
%V_l vorher schätzen
V_l_m=[41.967059865700620,  4.794498579122601, 1.589336549493195, 0.650927769700661,    0.292478109935620, 0.142717819629993, 0.069678551240208,  0.034946040052114, 0.016939479860485, 0.008511726065860, 0.004284090378739 ]; %Milstein V_l
V_l=[40.892827019950865, 5.121167149569525, 2.325199583324076, 0.986996681209547, 0.431678536189194, 0.198008343730391, 0.094348272635517, 0.046142300800308, 0.022590165062293, 0.011245767899094, 0.005705397165832, 0.002911580313002]; %Euler V_l
V_l_i=[41.369871034796262,1.574348141206355, 0.596726677655198, 0.255480890604358,0.116753536476469, 0.054953738733514, 0.027708045697653,0.013679815924426, 0.006768758270061, 0.003630122753334, 0.001577908995141]; %IJK V_l

L= optimale_L_asia_euler(V_l, epsilon, K, T, r, v0, kappa, theta, sigma, rho, M, M0, S0); %L für jedes Verfahren einzeln bestimmen
L_m= optimale_L_asia_milstein(V_l_m, epsilon, K, T, r, v0, kappa, theta, sigma, rho, M, M0, S0);
L_ijk=optimale_L_asia_ijk(V_l_i, epsilon, K, T, r, v0, kappa, theta, sigma, rho, M, M0, S0);

h_l = zeros(1, L+1);  
N_opt=zeros(L+1, length(epsilon));
sum_term=0;
    for l=0:L
         h_l(l+1) = T / (M0*M^l);
         sum_term = sum_term+ sqrt(V_l(l+1)/ (T/(M0*M^l)));
    end
    % Kosten pro Level proportional zu Zeitschritten; C_l=M^l
        for i = 1:length(epsilon)  
                eps=epsilon(i);
                for l=0:L
                N_opt(l+1, i) = ceil((1 / eps^2) * sqrt(V_l(l+1).* h_l(l+1)) * sum_term); % Vektor mit optimalen N_l % ceil aufrunden
                end
        end
h_l = zeros(1, L_m+1);  
N_opt_m=zeros(L_m+1, length(epsilon));
sum_term=0;
     for l=0:L_m
          h_l(l+1) = T / (M0*M^l);
          sum_term = sum_term+ sqrt(V_l_m(l+1)/ (T/(M0*M^l)));
     end
         for i = 1:length(epsilon)  
                 eps=epsilon(i);
                 for l=0:L_m
                      N_opt_m(l+1, i) = ceil((1 / eps^2) * sqrt(V_l_m(l+1).* h_l(l+1)) * sum_term); % Vektor mit optimalen N_l % ceil aufrunden
                 end
         end
h_l = zeros(1, L_ijk+1);  
N_opt_i=zeros(L_ijk+1, length(epsilon));
sum_term=0;
     for l=0:L_ijk
          h_l(l+1) = T / (M0*M^l);
          sum_term = sum_term+ sqrt(V_l_i(l+1)/ (T/(M0*M^l)));
     end
          for i = 1:length(epsilon)  
                  eps=epsilon(i);
                  for l=0:L_ijk
                       N_opt_i(l+1, i) = ceil((1 / eps^2) * sqrt(V_l_i(l+1).* h_l(l+1)) * sum_term); % Vektor mit optimalen N_l % ceil aufrunden
                  end
          end

figure('Units','normalized','Position',[0.1 0.1 0.7 0.6]);  % größere Figur
colors = lines(3);  % drei unterscheidbare Farben

% Euler
plot(0:L, N_opt, 'o-', ...
    'Color', colors(1,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(1,:), ...
    'MarkerEdgeColor', 'none');
hold on;

% Milstein
plot(0:L_m, N_opt_m, 'o-', ...
    'Color', colors(2,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(2,:), ...
    'MarkerEdgeColor', 'none');

% IJK
plot(0:L_ijk, N_opt_i, 'o-', ...
    'Color', colors(3,:), ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', colors(3,:), ...
    'MarkerEdgeColor', 'none');

xlabel('Level l');
ylabel(' Optimale N_l');
legend('Euler', 'Milstein', 'IJK')
title('Optimale Pfadanzahl pro Level für MLMC');
grid on;
end