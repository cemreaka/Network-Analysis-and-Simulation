%% Q2) 

clear all;
clc;

% a)

n = 48; 
rng('default');
rv = rand(1, n);

% b)

s_mean = mean(rv);
s_std = std(rv);
alpha = 0.05;
df = n - 1;
t_crit = tinv(1 - alpha / 2, df);
s_error = s_std / sqrt(n);
ci = s_mean + [-1 1] * t_crit * s_error;

fprintf('Sample Mean: %f\n', s_mean);
fprintf('Sample STD: %f\n', s_std);
fprintf('95%% Confidence Interval for Mean: [%f, %f]\n', ci(1), ci(2));

% c)

N = 1000;
res = zeros(N, 3);
for i = 1:N
    rv = rand(1, n);
    s_mean = mean(rv);
    s_std = std(rv);
    s_error = s_std / sqrt(n);
    ci = s_mean + [-1 1] * t_crit * s_error;
    true_mean = 0.5;
    contains_mean = (ci(1) <= true_mean) && (true_mean <= ci(2));
    res(i,:) = [s_mean, ci(1), ci(2)];
    if ~contains_mean
        res(i,3) = -1;
    end
end

sorted_res = sortrows(res, 2);
no_not_contain = sum(sorted_res(:, 3) == -1);
fprintf('No. of times the confidence interval does not contain the true value of the mean: %d out of %d\n', no_not_contain, N);

figure;
errorbar(sorted_res(:,1), sorted_res(:,2), sorted_res(:,3) - sorted_res(:,1), 'o');
xlabel('Sample mean');
ylabel('Confidence interval');
title('Sample means and confidence intervals');

%% Q4)

clear all;
clc;

% a)
n = 1000;
rv = rand(1, n);
s_mean = mean(rv);
s_var = var(rv);

fprintf('Sample mean: %f\n', s_mean);
fprintf('Sample variance: %f\n', s_var);

% b) General approach
n_max = 10000;
s_means = zeros(1, n_max);
s_vars = zeros(1, n_max);
mean_errors = zeros(1, n_max);
var_errors = zeros(1, n_max);

for i = 1:n_max
    rv = rand(1, i);
    s_means(i) = mean(rv);
    s_vars(i) = var(rv);
    mean_errors(i) = abs(s_means(i) - 0.5);
    var_errors(i) = abs(s_vars(i) - (1 / 12)); % (b - a)^2 / 12 
end

figure;
semilogx(1:n_max, mean_errors);
hold on;
semilogx(1:n_max, var_errors);
xlabel('n');
ylabel('Error');
title('Accuracy of sample mean and variance estimates');
legend('Mean error', 'Variance error');

% c) 
z = 1.96;
for i = 1:n_max
    ci(i,:) = [s_vars(i) - z * sqrt(s_vars(i) / i), s_vars(i) + z * sqrt(s_vars(i) / i)];
end

figure;
semilogx(1:n_max, var_errors);
hold on;
plot(1:n_max, ci(:,1), '--r');
plot(1:n_max, ci(:,2), '--r');
xlabel('Sample size (n)');
ylabel('Absolute error in variance estimate');
title('Accuracy of variance estimate for U(0,1) distribution');
legend('Error', '95% Confidence Interval');

% d)


%% Q5) For Q2

clear all;
clc;

% a)

n = 48; 
rng('default');
rv = randn(1, n);

% b)

s_mean = mean(rv);
s_std = std(rv);
alpha = 0.05;
df = n - 1;
t_crit = tinv(1 - alpha / 2, df);
s_error = s_std / sqrt(n);
ci = s_mean + [-1 1] * t_crit * s_error;

fprintf('Sample Mean: %f\n', s_mean);
fprintf('Sample STD: %f\n', s_std);
fprintf('95%% Confidence Interval for Mean: [%f, %f]\n', ci(1), ci(2));

% c)

N = 1000;
res = zeros(N, 3);
for i = 1:N
    rv = randn(1, n);
    s_mean = mean(rv);
    s_std = std(rv);
    s_error = s_std / sqrt(n);
    ci = s_mean + [-1 1] * t_crit * s_error;
    true_mean = 0.5;
    contains_mean = (ci(1) <= true_mean) && (true_mean <= ci(2));
    res(i,:) = [s_mean, ci(1), ci(2)];
    if ~contains_mean
        res(i,3) = -1;
    end
end

sorted_res = sortrows(res, 2);
no_not_contain = sum(sorted_res(:, 3) == -1);
fprintf('No. of times the confidence interval does not contain the true value of the mean: %d out of %d\n', no_not_contain, N);

figure;
errorbar(sorted_res(:,1), sorted_res(:,2), sorted_res(:,3) - sorted_res(:,1), 'o');
xlabel('Sample mean');
ylabel('Confidence interval');
title('Sample means and confidence intervals');

%% Q5) For Q4

clear all;
clc;

% a)
n = 1000;
rv = randn(1, n);
s_mean = mean(rv);
s_var = var(rv);

fprintf('Sample mean: %f\n', s_mean);
fprintf('Sample variance: %f\n', s_var);

% b) General approach
n_max = 10000;
s_means = zeros(1, n_max);
s_vars = zeros(1, n_max);
mean_errors = zeros(1, n_max);
var_errors = zeros(1, n_max);

for i = 1:n_max
    rv = randn(1, i);
    s_means(i) = mean(rv);
    s_vars(i) = var(rv);
    mean_errors(i) = abs(s_means(i) - 0.5);
    var_errors(i) = abs(s_vars(i) - (1 / 12)); % (b - a)^2 / 12 
end

figure;
semilogx(1:n_max, mean_errors);
hold on;
semilogx(1:n_max, var_errors);
xlabel('n');
ylabel('Error');
title('Accuracy of sample mean and variance estimates');
legend('Mean error', 'Variance error');

% c) 
z = 1.96;
for i = 1:n_max
    ci(i,:) = [s_vars(i) - z * sqrt(s_vars(i) / i), s_vars(i) + z * sqrt(s_vars(i) / i)];
end

figure;
semilogx(1:n_max, var_errors);
hold on;
plot(1:n_max, ci(:,1), '--r');
plot(1:n_max, ci(:,2), '--r');
xlabel('Sample size (n)');
ylabel('Absolute error in variance estimate');
title('Accuracy of variance estimate for U(0,1) distribution');
legend('Error', '95% Confidence Interval');

% d)