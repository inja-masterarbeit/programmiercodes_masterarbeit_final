function N_opt = optimale_Nl(V_l, M, epsilon, T, M0)

%M = 2;  % Verfeinerungsfaktor
%epsilon = [0.005, 0.01, 0.02,  0.05, 0.1];   % gewünschter RMSE
%T=1;
%V_l schätzen

L = length(V_l) - 1;           % Level von 0 bis L
N_opt = zeros(length(epsilon), L+1);
C_l = zeros(1, L+1);           % Kosten pro Level

% Kosten pro Level berechnen
C_l(1) = M0;                   % Level 0
for l = 2:L+1
    C_l(l) = M0*M^(l-1) + M0*M^(l-2); % Level l>=1
end

sum_term = sum(sqrt(V_l .* C_l));

% N_l für jede gewünschte Fehlertoleranz
for i = 1:length(epsilon)
    eps = epsilon(i);
    for l = 0:L
        N_opt(i, l+1) = ceil( (1 / eps^2) * sqrt(V_l(l+1)/C_l(l+1)) * sum_term );
    end
end

end
