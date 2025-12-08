clear
clc
%% Load data (Matlab Grader provides getSPOdata)
birthdate = 20030509;               % yyyymmdd, given in uppgiften
[lambdavec, Tvec, cvec] = getSPOdata(birthdate);

%% Question 1:
all_EBO = zeros(1,11);
for j=1:9 %All LRUs 

    mu = lambdavec(j) * Tvec(j);
    
    for s=0:11 % All amount of spare parts
        all_EBO(s+1) = EBO(s, mu);
    end

end

if all(diff(all_EBO))
    disp("It's integer-convex")
end


%% Question 2:
n = length(lambdavec);  
EBO0_j = zeros(size(lambdavec));
EBO1_j = zeros(size(lambdavec));

for j = 1:n
    mu = lambdavec(j) * Tvec(j);
    EBO0_j(j) = EBO(0, mu);
    EBO1_j(j) = EBO(1, mu);
end

EBO0_total = sum(EBO0_j);
EBO2_candidates = EBO0_total - EBO0_j + EBO1_j;

[EBO2, idx] = min(EBO2_candidates);
Cost2 = cvec(idx);
Q2 = [EBO2, Cost2];


%% Question 3:
n = length(lambdavec);
Cmax = 500;
mu = lambdavec .* Tvec;

% Start with no spare parts
s_curr = zeros(1,n);
Cost_curr = 0;
EBO_curr = f(s_curr, mu);

S_history = s_curr;
Cost_history = Cost_curr;
EBO_history = EBO_curr;

% MALLOC Algorithm
while true
    best_ratio = 0;
    best_j = 0;
    best_EBO_new = EBO_curr;

    for j=1:n
        if Cost_curr + cvec(j) <= Cmax
            s_try = s_curr;
            s_try(j) = s_try(j) +1;

            EBO_new = f(s_try, mu);
            deltaEBO = EBO_curr - EBO_new;

            if deltaEBO > 0
                ratio = deltaEBO / cvec(j);
                if ratio > best_ratio
                    best_ratio = ratio;
                    best_j = j;
                    best_EBO_new = EBO_new;
                end
            end
        end
    end
    if best_j == 0
        break;
    end

    s_curr(best_j) = s_curr(best_j) +1;
    Cost_curr = Cost_curr + cvec(best_j);
    EBO_curr = best_EBO_new;

    S_history = [S_history; s_curr];
    Cost_history = [Cost_history; Cost_curr];
    EBO_history = [EBO_history; EBO_curr];
end


%% Uppgift 4
figure;
plot(Cost_history, EBO_history, '-o')
xlabel('Total Cost')
ylabel('Expected Backorders');
title('Efficient Solutions Curve - Marginal Allocation')
EPtable = [ S_history(1:5, :) , EBO_history(1:5), Cost_history(1:5) ];

%% Uppgift 6
% I rapporten

%% Uppgift 7
Cmax1 = 50;
j = 1;                       % LRU1
mu1 = lambdavec(j) * Tvec(j);
c1  = cvec(j);


smax1 = floor(Cmax1 / c1);

EBO1 = zeros(smax1+1,1);
for s = 0:smax1                 % loopar upp till 3
    EBO1(s+1) = EBO(s, mu1);
end

bestEBO = zeros(Cmax1+1,1);
bestS   = zeros(Cmax1+1,1);   
for b = 0:Cmax1
    best_val = inf;
    best_s   = 0;

    max_s_for_b = floor(b / c1);
    for s = 0:max_s_for_b
        val = EBO1(s+1); 
        if val < best_val
            best_val = val;
            best_s   = s;
        end
    end

    bestEBO(b+1) = best_val;
    bestS(b+1)   = best_s;
end

LRU1 = bestS.';   % 1×51
%% Uppgift 8

B = 500;
n = length(lambdavec);
mu = lambdavec(:) .* Tvec(:);
c = cvec(:);

smax = floor(B ./ c);
EBOtable = cell(n,1);
for j = 1:n
    EBO_j = zeros(smax(j)+1, 1);
    for s = 0:smax(j)
        EBO_j(s+1) = EBO(s, mu(j));
    end
    EBOtable{j} = EBO_j;
end

% DP-tabellen V (n+1) x (B+1)
% DP-tabellen V (n+1) x (B+1)
% V(k, b+1) = min EBO för LRU k..n med budget b
V = zeros(n+1, B+1);
policy = zeros(n, B+1);    % policy(k, b+1) = optimalt s_k vid budget b

% terminal: V(n+1, :) = 0 (redan nollor)

for k = n:-1:1                      % från sista LRU till första
    smax_k = smax(k);
    ck     = c(k);
    EBO_k  = EBOtable{k};           % vektor EBO_k(0..smax_k)

    for b = 0:B
        best_val = inf;
        best_s   = 0;
        max_s_for_b = min(smax_k, floor(b / ck));

        for s = 0:max_s_for_b
            val = EBO_k(s+1) + V(k+1, b - ck*s + 1);  

            if val < best_val
                best_val = val;
                best_s   = s;
            end
        end

        V(k, b+1)      = best_val;
        policy(k, b+1) = best_s;
    end
end

DP_table = zeros(B+1, n);  
EBO_DP   = zeros(B+1, 1);
Cost_DP  = zeros(B+1, 1);

for B_now = 0:B
    b = B_now;
    s = zeros(1,n);

    for k = 1:n
        s_k = policy(k, b+1);
        s(k) = s_k;
        b = b - c(k)*s_k;
    end

    DP_table(B_now+1, :) = s;

    totalEBO = 0;
    for j = 1:n
        totalEBO = totalEBO + EBO(s(j), mu(j));
    end
    EBO_DP(B_now+1)  = totalEBO;
    Cost_DP(B_now+1) = s * cvec(:); 
end

budgets = [0 100 150 350 500];
numB    = length(budgets);
DynPtable = zeros(numB, n+2);

for i = 1:numB
    B_now = budgets(i);
    DynPtable(i, :) = [DP_table(B_now+1, :), EBO_DP(B_now+1), Cost_DP(B_now+1)];
end


% Hämta lösningar
budgets = [0 100 150 350 500];
numB    = length(budgets);
X_opt   = zeros(numB, n);
EBO_opt = zeros(numB, 1);
Cost_opt= zeros(numB, 1);

for idx = 1:numB
    B_now = budgets(idx);
    b = B_now;
    s = zeros(1,n);

    for k = 1:n
        s_k = policy(k, b+1);
        s(k) = s_k;
        b = b - c(k)*s_k;
    end

    X_opt(idx,:) = s;

    totalEBO = 0;
    for j = 1:n
        totalEBO = totalEBO + EBO(s(j), mu(j));
    end
    EBO_opt(idx)  = totalEBO;
    Cost_opt(idx) = X_opt(idx,:)*c;
end

DynPtable = [X_opt, EBO_opt, Cost_opt];

%% Uppgift 9
b_vec = 0:B;
EBO_DP = V(1, :);

figure;
plot(b_vec, EBO_DP, 'LineWidth', 1.5); 
hold on;

plot(Cost_history, EBO_history, 'LineWidth', 1.5);

hold off;
grid on;
xlabel('Total cost / budget');
ylabel('Expected backorders (EBO)');
title('Dynamic Programming vs Marginal Allocation');
legend('Dynamic Programming (DP)', 'Marginal Allocation', 'Location', 'northeast');

%% Uppgift 10 – jämför mot "en av varje LRU"
s_one   = ones(1,n);
EBO_one = f(s_one, lambdavec .* Tvec);
Cost_one = s_one * cvec(:);


mask_a = (EBO_DP <= EBO_one + 1e-8);
idx_a  = find(mask_a);

[Cost_best_a, pos_a] = min(Cost_DP(idx_a)); 
best_a = idx_a(pos_a);                      

x_a   = DP_table(best_a, :);                
EBO_a = EBO_DP(best_a);                     

Costimprovement = Cost_one - Cost_best_a;  


mask_b = (Cost_DP <= Cost_one + 1e-8);   
idx_b  = find(mask_b);

[EBO_best_b, pos_b] = min(EBO_DP(idx_b)); 
best_b = idx_b(pos_b);

x_b    = DP_table(best_b, :);               
Cost_b = Cost_DP(best_b);                

EBOimprovement = EBO_one - EBO_best_b; 

Q10 = [EBOimprovement, Costimprovement]
%% Functions  
function exp_backord = EBO(s, mu)
    k_max = s + 30 + round(10*mu);
    ks = (s+1):k_max;
    probs = poisspdf((s+1):k_max, mu);
    exp_backord = sum((ks - s) .* probs);
end

function totalEBO = f(s, muvec)
    % totalEBO(s, muvec) = sum_j EBO(s_j, mu_j)
    totalEBO = 0;
    for j = 1:length(s)
        totalEBO = totalEBO + EBO(s(j), muvec(j));
    end
end