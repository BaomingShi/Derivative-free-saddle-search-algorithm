rng(1)
d=2;
xstar=zeros(d,1);

E=@(x) slightly_harder_saddle(x);

u=[0.8;0.8];
options=struct('k',1,'innpfunc',@(x,y)(x'*y),'maxiter',5000,'dt',0.001);options.outputX=1;options.gammamax=10;options.gammamin=0.5;options.betat=5;options.betau=0.2;
options.V=randn(d,options.k);
options.V=options.V(:,1:options.k);
options.stepsize1=0.01;
options.l=0.1;


options.decaystep=0;
options.stepsize2=100;


tic
[u, x_total,output] = stochastic_hiosd_sirqit_zerothorder(E,u,options);
t=toc;
fprintf('The CPU time of zeroth-order algorithm is %f seconds\n', t);


[d,numberiter]=size(x_total);
error=zeros(1,numberiter);
for i=1:numberiter
    error(i)=norm(x_total(:,i)-xstar).^2;
end

figure(1)
[xdata,ydata]=meshgrid([-1:0.01:1],[-1:0.01:1]);

for i=1:length([-1:0.01:1])
    for j=1:length([-1:0.01:1])
        zdata(i,j)=E([xdata(i,j);ydata(i,j)]);
    end
end
contour(xdata,ydata,zdata,50,'HandleVisibility','off');
hold on
plot(x_total(1,1:1:end),x_total(2,1:1:end),'--*','color',[1 0.1 0.3],'LineWidth',1)

scatter(0, 0, 50, [0, 0, 0], 'o', 'filled', 'LineWidth', 2);

scatter(0.8, 0.8, 50, [0.1, 0.6, 0.8], 'o', 'filled', 'LineWidth', 2);

ll=legend('Trajectory','Saddle point $\mathbf{x}^*$', 'Initial condition $\mathbf{x}(0)$','inter');
set(ll,'interpreter','latex')
set(gca,'fontsize',16)

figure(2)
semilogy(1:1:numberiter,error,'color','black','LineWidth',1)
hold on

xlabel('n')
ylabel('$\Vert \mathbf{x}(n)-\mathbf{x}^* \Vert^2$','Interpreter','latex')
set(gca,'fontsize',16)
grid on



function val = slightly_harder_saddle(x)
    lr = 0.2;      
    steps = 30;     
    z=[0.1;0.05];
    kk=2;
    for i = 1:steps
        grad = [kk*(z(1)-x(1))^(kk-1) + z(2)*cos(z(1)*z(2));
                kk*(z(2)-x(2))^(kk-1) + z(1)*cos(z(1)*z(2))];  % ∇_z g
        z = z - lr * grad;
    end
    
    val = (z(1)-x(1))^kk + (z(2)-x(2))^kk + sin(z(1)*z(2));
end

