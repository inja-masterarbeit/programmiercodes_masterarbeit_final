function costs=gesamtkosten_vs_epsilon (V_l, V_l_m, V_l_i, epsilon,M,T,M0)
M=2;M0=4;T=1;
epsilon = [0.0001, 0.0005, 0.005, 0.01, 0.02,  0.05]; % gewünschter RMSE

V_l_i=[41.369871034796262,1.574348141206355, 0.596726677655198, 0.255480890604358,0.116753536476469, 0.054953738733514, 0.027708045697653,0.013679815924426, 0.006768758270061, 0.003630122753334, 0.001577908995141]; %IJK V_l
V_l=[40.892827019950865, 5.121167149569525, 2.325199583324076, 0.986996681209547, 0.431678536189194, 0.198008343730391, 0.094348272635517, 0.046142300800308, 0.022590165062293, 0.011245767899094, 0.005705397165832, 0.002911580313002]; %Euler V_l
V_l_m=[41.967059865700620,  4.794498579122601, 1.589336549493195, 0.650927769700661,    0.292478109935620, 0.142717819629993, 0.069678551240208,  0.034946040052114, 0.016939479860485, 0.008511726065860, 0.004284090378739 ]; %Milstein V_l
V_l_i_cv=[ 3.452893435625576, 0.692929053896580, 0.375500580076401,0.198512331038281,  0.104177619199234, 0.052366843968889, 0.026402115418500, 0.013380749250153, 0.006564654703669, 0.003289595957638,  0.001656114673772, 0.000830315331419];
V_l_e_cv=[  7.025755516530939,  4.615852617098510,  2.030390330691740, 0.812464081598686, 0.327687296408386, 0.144685323935348, 0.067627253894367, 0.032605409720654, 0.015956735664344, 0.007899452462152,  0.004000387051057];
V_l_m_cv=[6.907216688468123, 4.293824930604254,  1.449561395288868, 0.617554223292023, 0.298001097221532, 0.143654946703399, 0.071967081664341,    0.036035828687824,0.017936915858768, 0.008793845731158, 0.004539719888657, 0.002288614543551];

N_l = optimale_Nl(V_l, M, epsilon, T, M0);
N_l_i = optimale_Nl(V_l_i, M, epsilon, T, M0);
N_l_m = optimale_Nl(V_l_m, M, epsilon, T, M0);
N_l_i_cv = optimale_Nl(V_l_i_cv, M, epsilon, T, M0);
N_l_e_cv = optimale_Nl(V_l_e_cv, M, epsilon, T, M0);
N_l_m_cv = optimale_Nl(V_l_m_cv, M, epsilon, T, M0);

C_l=zeros(length(epsilon),length(V_l)); %Kosten
C_l_i=zeros(length(epsilon),length(V_l));
C_l_m=zeros(length(epsilon),length(V_l));
C_l_i_cv=zeros(length(epsilon),length(V_l));
C_l_e_cv=zeros(length(epsilon),length(V_l));
C_l_m_cv=zeros(length(epsilon),length(V_l));

L_e=[9,8,8,7,6,4]; %optimalen L vorab berechnet (auch über optimale_L Code möglich)

for i=1:length(epsilon)
    L=L_e(i);

    for l=0:L
 
       C_l(i,l+1)=N_l(i,l+1)*M0*M^l; %Kosten auf dem Level l, Euler
       C_l_i(i,l+1)=N_l_i(i,l+1)*M0*M^l;
       C_l_m(i,l+1)=N_l_m(i,l+1)*M0*M^l;
       C_l_i_cv(i,l+1)=N_l_i_cv(i,l+1)*M0*M^l;
       C_l_e_cv(i,l+1)=N_l_e_cv(i,l+1)*M0*M^l;
       C_l_m_cv(i,l+1)=N_l_m_cv(i,l+1)*M0*M^l;
  
    end
    
end
    costs=log(sum(C_l,2));
    costs_i=log(sum(C_l_i,2));
    costs_m=log(sum(C_l_m,2));
    costs_i_cv=log(sum(C_l_i_cv,2));
    costs_e_cv=log(sum(C_l_e_cv,2));
    costs_m_cv=log(sum(C_l_m_cv,2));
    disp(costs_i);
    disp(costs);

    x=log(epsilon);
    y=log(epsilon.^(-2).*log(epsilon).^2*500); %Referenzlinie für Fall beta=1, ohne cv
   
    y1=log(epsilon.^(-3)*500); 

figure; hold on;
hold on
h1 = plot(costs,   x, '--d');
h2 = plot(costs_m, x, '--d');
h3 = plot(costs_i, x, '--d');
h4 = plot(y,  x, 's:');
h5 = plot(y1, x, 's:');
h6 = plot(costs_e_cv, x, '-o', 'LineWidth', 2);
h7 = plot(costs_m_cv, x, '-o', 'LineWidth', 2);
h8 = plot(costs_i_cv, x, '-o', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);

% Alle Linien in einem Array
H = [h1 h2 h3 h4 h5 h6 h7 h8];

% Markerfarbe = Linienfarbe, keine Kontur
for i = 1:numel(H)
    c = get(H(i), 'Color');                     % Linienfarbe auslesen
    set(H(i), 'MarkerFaceColor', c, ...         % Füllung = Linienfarbe
               'MarkerEdgeColor', 'none');      % keine Kontur
end
colors = lines(8);  % oder: parula(8), turbo(8), hsv(8), etc.
H = [h1 h2 h3 h4 h5 h6 h7 h8];

for i = 1:numel(H)
    set(H(i), 'Color', colors(i,:), ...
               'MarkerFaceColor', colors(i,:), ...
               'MarkerEdgeColor', 'none');
end
hold off
hold off;
p = polyfit(log(epsilon), costs, 1);
rate = -p(1);
disp(rate);
ylabel('$\varepsilon$ ', 'Interpreter', 'latex','FontWeight','bold','FontSize',14);
xlabel('Gesamtkosten C ');


title('Gesamtkosten der Verfahren in Abhängigkeit der Genauigkeit ε (log-log)', 'Interpreter', 'latex');
legend({'Euler', 'Milstein', 'IJK', '$\varepsilon^{-2} \cdot \log(\varepsilon)^2$', ...
        '$\varepsilon^{-3}$', 'Euler+CV', 'Milstein+CV', 'IJK+CV'}, ...
        'Interpreter','latex');
grid on;
 
 

end

