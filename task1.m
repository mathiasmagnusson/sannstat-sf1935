clc; clear variables; close all;

alpha = 2;
sigmasq = 0.6;
beta = 1 / sigmasq;
[w0, w1] = meshgrid(-2:0.01:2);
W = [w0(:), w1(:)];

prior = reshape(mvnpdf(W, 0, 1/alpha * eye(2)), length(w1), length(w0));
set(0, defaultTextInterpreter="latex", defaultLegendInterpreter="latex", defaultAxesTickLabelInterpreter="latex");
figure(Name="Prior distribution")
contourf(w0, w1, prior, 10)
xlabel("$w_0$", FontSize=15);
ylabel("$w_1$", FontSize=15);
axis equal

w = [-1.5; 0.5];
x = -1:0.01:1;
x_ext = [ones(size(x)); x];
%rng("default");
epsilon = normrnd(0, sqrt(sigmasq), size(x));
t = w' * x_ext + epsilon;
figure(Name="Data set")
scatter(x, t, 25, "filled")
xlabel("$x$", FontSize=15);
ylabel("$t$", FontSize=15);


% for i = 1
%     likelihood = reshape(normpdf(t(i), W * x_ext(:, i), 1 / beta), length(w1), length(w0));
%     posterior = likelihood .* prior;
%     prior = posterior;
% end

figure;
tiledlayout(7,3)
subset = randperm(numel(t));
idx = subset(1:7);
for i = 1:numel(idx)
    likelihood = reshape(normpdf(t(idx(i)), W * x_ext(:, idx(i)), 1 / beta), length(w1), length(w0));
    posterior = likelihood .* prior;
    prior = posterior;
%     if i < 99
%         continue
%     end

    nexttile
    contourf(w0, w1, likelihood, 10)
    %title("Likelihood", FontSize=16)
    xlabel("$w_0$", FontSize=15);
    ylabel("$w_1$", FontSize=15);
    
    
    
    nexttile;
    contourf(w0, w1, posterior, 10)
    %title("Posterior", FontSize=16)
    xlabel("$w_0$", FontSize=15);
    ylabel("$w_1$", FontSize=15);
    
    


    nexttile
    Sn = inv(alpha * eye(2) + beta * (x_ext(:, idx(1:i)) * x_ext(:, idx(1:i))'));
    mn = beta * Sn * x_ext(:, idx(1:i)) * t(:, idx(1:i))';
    plot(x, mvnrnd(mn, Sn, 5) * x_ext)
    hold on;
    scatter(x, t, 25, "filled")
    scatter(x(idx(1:i)), t(idx(1:i)), 25, "k", "filled")
    hold off;
end



% %nexttile
% %figure(Name="Prior distribution")
% contourf(w0, w1, prior, 10)
% xlabel("$w_0$", FontSize=15);
% ylabel("$w_1$", FontSize=15);
% title("Prior distribution", FontSize=16)
% axis equal
% 
% figure(Name="Likelihood")
% %nexttile
% contourf(w0, w1, likelihood, 10)
% xlabel("$w_0$", FontSize=15);
% ylabel("$w_1$", FontSize=15);
% title("Likelihood", FontSize=16)
% axis equal
% %nexttile
% figure(Name="Posterior distribution")
% contourf(w0, w1, posterior, 10)
% xlabel("$w_0$", FontSize=15);
% ylabel("$w_1$", FontSize=15);
% title("Posterior distribution", FontSize=16)
% axis equal
% 
% subset = randperm(numel(t));
% idx = subset(1:7);
% 
% tiledlayout("flow")
% 
% for i = 1:numel(idx)
%     nexttile
%     Sn = inv(alpha * eye(2) + beta * (x_ext(:, idx(1:i)) * x_ext(:, idx(1:i))'));
%     mn = beta * Sn * x_ext(:, idx(1:i)) * t(:, idx(1:i))';
%     plot(x, mvnrnd(mn, Sn, 5) * x_ext)
%     hold on;
%     title("$n=" + i + "$", FontSize=12)
%     plot(x, t, ".")
% end


