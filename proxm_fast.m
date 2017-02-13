function x = proxm_fast( y,m )
%   Prox onto x'1 = m and x_k \in [0,1] \forall k
n = length(y);
[v,i] = sort(y,'ascend');

% first figure out how many are going to be 0
val_m = v * ones(1,n) - ones(n,1) * v';
val_m_crop = max(min(val_m,1),0);

val_sum = sum(val_m_crop);
k0 = sum(val_sum >= m);

% then figure out how many are going to be 1
y1 = -y;
[v1,~] = sort(y1,'ascend');
val_m1 = (v1-1) * ones(1,n) - ones(n,1) * v1';
val_m1_crop = max(min(val_m1,0),-1);

val_sum_1 = sum(val_m1_crop);
k1 = sum(val_sum_1 > -m);


% now we can solve for lambda
n_between = n - k0 - k1;
if n_between > 0
    lambda = (m - k1 - sum(v(k0+1:n-k1))) / n_between;
    xmag = max(min(v + lambda,1),0);
else
    xmag = [zeros(k0,1); ones(k1,1)];
end


x = zeros(n,1);
x(i) = xmag;

end

