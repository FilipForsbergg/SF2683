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
EPtable = [ S_history(1:5, :) , EBO_history(1:5), Cost_history(1:5) ]


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