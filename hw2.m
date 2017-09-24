alpha=0.35
beta=0.99
delta=0.025
sigma=2

pi=[0.977,0.023; 0.926,0.074]
A=[1.1;0.678]
gridpoints=1000
shocks=2
convergence=1
tolerance=1e-06

k_min = 0;
k_max = 1.1*(alpha*A(2)/(1/beta-1+delta))^(alpha/(1-alpha))+(1-delta)*(alpha*1.1/(1/beta-1+delta))^(1/(1-alpha));
k=linspace(k_min,k_max,gridpoints)
v=zeros(2, gridpoints)
v_guess = zeros(2, gridpoints);
consumption=zeros(2, gridpoints)
utility=zeros(2,gridpoints)

while convergence>tolerance
    for i=1:gridpoints
        for j=1:shocks
            consumption=A(j)*k(i)^alpha+(1-delta)*k(i)-k
            minus=find(consumption<0)
            consumption(minus)=NaN
            utility(j,:)=(consumption.^(1-sigma))/(1-sigma)
            utility(j,minus)=-1e12
        end
        [v_guess(:,i),p(:,1)]=max(utility+beta*(v*pi))
    end
    con=max(abs(v_guess-v))
    v_guess=v
    i=i+1
end

plot(v(1,:))
figure(plot(v(2,:)))            