
rng(1)
d=2;
xstar=ones(d,1);
s=ones(d,1);
s(1:1)=-50;

[~, ~, EValue] = get_Rm_handles(s);

E=@(x) (sum(100*(x(2:end)-x(1:end-1).^2).^2+(1-x(1:end-1)).^2)+sum(s.*atan(x-xstar).^2));


u=xstar-0.5*randn(d,1);
options=struct('k',1,'innpfunc',@(x,y)(x'*y),'maxiter',20000,'dt',0.001);options.outputX=1;options.gammamax=10;options.gammamin=0.5;options.betat=5;options.betau=0.2;

options.V=randn(d,1);
options.l=0.0001;

options.stepsize1=0.00001;


options.decaystep=0;
options.stepsize2=100;


tic
[u_end, x_total,output] = stochastic_hiosd_sirqit_zerothorder(E,u,options);
t=toc;
fprintf('The CPU time of zeroth-order algorithm is %f seconds\n', t);


[d,numberiter]=size(x_total);
error=zeros(1,numberiter);
for i=1:numberiter
    error(i)=norm(x_total(:,i)-xstar).^2;
end

figure(1)
semilogy(1:1:numberiter,error,'color','black','LineWidth',1)
hold on

xlabel('n')
ylabel('$\Vert \mathbf{x}(n)-\mathbf{x}^* \Vert^2$','Interpreter','latex')
set(gca,'fontsize',16)
grid on


figure(2)
[x1, x2] = meshgrid(linspace(0, 2, 200), linspace(0, 2, 200));
z = zeros(size(x1));


for i = 1:numel(x1)
    z(i) = E([x1(i);x2(i)]);
end



surf(x1, x2, z, 'EdgeColor', 'none','HandleVisibility','off');
colormap(pink); 
hold on;


x_points = x_total(1,:);
y_points = x_total(2,:);
z_points = arrayfun(@(a,b) E([a;b]), x_points, y_points, 'UniformOutput', false);
z_points = cell2mat(z_points);

plot3(x_points, y_points, z_points, 'ro-', ...
    'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'r');


xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
zlabel('$f(x_1,x_2)$','Interpreter','latex');
scatter3(1, 1, E([1;1])+10, 50, [0.8, 0.3, 0.1], 'o', 'filled', 'LineWidth', 2);

scatter3(u(1), u(2),E(u)+10 ,50, [0.1, 0.6, 0.8], 'o', 'filled', 'LineWidth', 2);

%ll=legend('Trajectory','Saddle point $\mathbf{x}^*$', 'Initial condition $\mathbf{x}(0)$','inter');
% set(ll,'interpreter','latex')
% set(gca,'fontsize',16)

grid on;
% view(45, 30);
% colorbar;


function [grad_handle, hess_handle, Rm_handle] = get_Rm_handles(s)
    Rm_handle = @(x) Rm_value(x, s);
    

    grad_handle = @(x) Rm_gradient(x, s);

    hess_handle = @(x) Rm_hessian(x, s);
end


function val = Rm_value(x, s)
    d = length(x);
    x_star = ones(d,1);
    R = 0;
    for i = 1:d-1
        R = R + 100*(x(i+1) - x(i)^2)^2 + (1 - x(i))^2;
    end
    arctan_part = atan(x - x_star);
    val = R + sum(s .* (arctan_part.^2));
end

% ----------------------------------------
function g = Rm_gradient(x, s)
    d = length(x);
    x_star = ones(d,1);
    g = zeros(d,1);
    

    for i = 1:d-1
        g(i) = g(i) - 400*(x(i+1) - x(i)^2)*x(i) - 2*(1 - x(i));
        g(i+1) = g(i+1) + 200*(x(i+1) - x(i)^2);
    end


    diff = x - x_star;
    arctan_term = atan(diff);
    g = g + 2 .* s .* arctan_term ./ (1 + diff.^2);
    g=-g;
end

% ----------------------------------------
function H = Rm_hessian(x, s)
    d = length(x);
    x_star = ones(d,1);
    H = zeros(d,d);
    

    for i = 1:d-1
        H(i,i) = H(i,i) + (-400*x(i+1) + 1200*x(i)^2 + 2);
        H(i,i+1) = H(i,i+1) - 400*x(i);
        H(i+1,i) = H(i+1,i) - 400*x(i);
        H(i+1,i+1) = H(i+1,i+1) + 200;
    end


    diff = x - x_star;
    H_arctan = diag(2 .* s .* (1 - diff.^2) ./ (1 + diff.^2).^2);
    
    H = H + H_arctan;
end
