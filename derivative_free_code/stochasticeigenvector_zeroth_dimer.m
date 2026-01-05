function v=stochasticeigenvector_zeroth_dimer(E,x,v)
k=size(v);
k=k(2);
Maxstep=10;

iteration=1;
F= @(uuu,d) (E(uuu+0.001*d)-E(uuu-0.001*d))/0.002*d;


while iteration<=Maxstep
    d=randn(size(x));
    Hv=(F(x+0.001*v(:,1),d)-F(x-0.001*v(:,1),d))/0.002;
    v(:,1)=v(:,1)-0.0002/length(x)*(Hv-Hv'*v(:,1)*v(:,1));
    v(:,1)=v(:,1)/norm(v(:,1));
    iteration=iteration+1;
end


for i=2:k
    iteration=1;
    v(:,i)=v(:,i)-v(:,1:i-1)*v(:,1:i-1)'*v(:,i);
    v(:,i)=v(:,i)/norm(v(:,i));
    while iteration<=Maxstep
        d=randn(size(x));
        Hv=(F(x+0.001*v(:,i),d)-F(x-0.001*v(:,i),d))/0.002;
        v(:,i)=v(:,i)-0.0002/length(x)*(Hv-v(:,1:i)*(v(:,1:i)'*Hv));
        v(:,i)=v(:,i)/norm(v(:,i));
        iteration=iteration+1;
    end
end


end