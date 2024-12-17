%% DDE Simulation
% Select Parameters
p.alpha = 1; p.beta = 60;
p.a = -1; p.b = -0.4;
p.c = -1; p.d = 0;
% p.tau1 = 0.076538085796; p.tau2 = p.tau1;
p.tau1 = 0.129; p.tau2 = p.tau1;

% Find u* and v* such that θu = 0.2 & θv = 0.2
p.theta_u = 0.7;
p.theta_v = 0.5;

biases = @(vars) [
        f_inv(vars(1),p) - (p.a .* vars(1)) - (p.b .* vars(2)) - p.theta_u;
        f_inv(vars(2),p) - (p.c .* vars(1)) - (p.d .* vars(2)) - p.theta_v;
    ];

sol = fsolve(biases, [0.5, 0.5], optimoptions('fsolve', 'Display', 'off'));

p.u = sol(1);
p.v = sol(2);

clear sol

% Define DDE parameters
p.tspan = [0 10];
p.delays = [p.tau1 p.tau2];
p.history = [0.074 0.077];
p.options = ddeset('RelTol', 1e-5);

% Solve DDE
sol = dde23(@(t, y, Z) ddefun(t, y, Z, p), ...
    p.delays,p.history,p.tspan,p.options);

% Initialise figure 4
figure;
clf; hold on

% Plot chotic solution
plot(sol.y(1,1), sol.y(2,1), 'k.', 'markersize', 12)
for i = 1:length(sol.y)-1
    plot(sol.y(1,i:i+1), sol.y(2,i:i+1), ...
        'k', 'linewidth', 1.5); % u against v
    % drawnow;
end

% Format axes
xlabel("$\mathit{u}$", 'Interpreter', 'latex')
ylabel("$\mathit{v}$", 'Interpreter', 'latex','rotation',0)
set(gca,'FontSize', 14, 'FontName', 'Times')

%% --------------------------------------------------------------------- %%
% ------------------------------- f(x,p) -------------------------------- %
    % Define the inverse sigmoid function, for z = u,v
    function f = f(z,p)
        f = 1 ./ (1 + exp(-p.beta * z));
    end

% ----------------------------- f_inv(x,p) ------------------------------ %
    % Define the inverse sigmoid function, for z = u,v
    function f = f_inv(z,p)
        f = (1 ./ p.beta) .* log(z ./ (1 - z));
    end

% ---------------------------- ddefun(t,y,Z) ---------------------------- %
    function d = ddefun(~,y,Z,p)
        dudt = -y(1) + ...
            f(p.theta_u + p.a .* Z(1,1) + ...
            p.b .* Z(2,2),p);
        dvdt = p.alpha .* ...
            (-y(2) + f(p.theta_v + p.c .* Z(1,2) + ...
            p.d .* Z(2,1),p));
    
        d = [dudt; dvdt];
    end