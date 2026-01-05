clc
clear
A=[-200, -100, -170, 15];
a = [-1, -1, -6.5, 0.7];
b = [0, 0, 11, 0.6];
c =[-10, -10, -6.5, 0.7];
x= [1, 0, -0.5, -1];
y = [0, 0.5, 1.5, 1];


E= @(uuu) sum(A.*exp(a.*(ones(1,4)*uuu(1)-x).^2+b.*(ones(1,4)*uuu(1)-x).*(ones(1,4)*uuu(2)-y)+c.*(ones(1,4)*uuu(2)-y).^2));
F= @(uuu,d) -(E(uuu+0.00001*d)-E(uuu-0.00001*d))/0.00002*d;



u=[0;1];


ustar=[-0.822001557452550;0.624312804305146];
options=struct('k',1,'innpfunc',@(x,y)(x'*y),'maxiter',1000,'dt',0.001);options.outputX=1;options.gammamax=10;options.gammamin=0.5;options.betat=5;options.betau=0.2;
options.V=u;
options.stepsize1=0.0001;
options.decaystep=0;
tic
options.l=0.001;
[u, x_total,output] = stochastic_hiosd_sirqit_zerothorder(E,u,options);
t=toc;
fprintf('The CPU time of zeroth-order algorithm is %f seconds\n', t);


F=@(uuu) -[sum(A.*exp(a.*(ones(1,4)*uuu(1)-x).^2+b.*(ones(1,4)*uuu(1)-x).*(ones(1,4)*uuu(2)-y)+c.*(ones(1,4)*uuu(2)-y).^2).*(2*a.*(ones(1,4)*uuu(1)-x)+b.*(ones(1,4)*uuu(2)-y)));  sum(A.*exp(a.*(ones(1,4)*uuu(1)-x).^2+b.*(ones(1,4)*uuu(1)-x).*(ones(1,4)*uuu(2)-y)+c.*(ones(1,4)*uuu(2)-y).^2).*(b.*(ones(1,4)*uuu(1)-x)+2*c.*(ones(1,4)*uuu(2)-y)))];   
E=@(uuu) sum(A.*exp(a.*(ones(1,4)*uuu(1)-x).^2+b.*(ones(1,4)*uuu(1)-x).*(ones(1,4)*uuu(2)-y)+c.*(ones(1,4)*uuu(2)-y).^2));



figure(1)
[xdata,ydata]=meshgrid([-1.4:0.01:0.2],[0:0.01:1.8]);

for i=1:length([0:0.01:1.8])
    for j=1:length([-1.4:0.01:0.2])
        zdata(i,j)=E([xdata(i,j);ydata(i,j)]);
    end
end
contour(xdata,ydata,zdata,50,'HandleVisibility','off');
hold on
plot(x_total(1,1:1:end),x_total(2,1:1:end),'--*','color',[1 0.1 0.3],'LineWidth',1)

scatter(-0.822001557452550, 0.624312804305146, 50, [0, 0, 0], 'o', 'filled', 'LineWidth', 2);

scatter(0, 1, 50, [0.1, 0.6, 0.8], 'o', 'filled', 'LineWidth', 2);

ll=legend('Trajectory','Saddle point $\mathbf{x}^*$', 'Initial condition $\mathbf{x}(0)$','inter');
set(ll,'interpreter','latex')
set(gca,'fontsize',16)

figure(2)
x_total=(x_total(1,:)-ustar(1)).^2+(x_total(2,:)-ustar(2)).^2;

M=length(x_total);
semilogy(1:1:M,x_total,'color','black','LineWidth',1)
hold on

xlabel('n')
ylabel('$\Vert \mathbf{x}(n)-\mathbf{x}^* \Vert^2$','Interpreter','latex')
set(gca,'fontsize',16)
grid on
