function norm_w = compute_nHbnorm(g,a, L,t0)
dw = @(x) compute_idu(a, x);
g_nHb = @(x) (L-t0).*g((L-t0)*x+t0);
w_func = @(x) (g_nHb(x) - dw(x)).^2;
norm_w = sqrt(sum(vect_int_sympson(w_func, 0, 1, 10000))); %% Norm of (H^{1,0}_B)^n
norm_w = intval(sup(norm_w));
end

function idu = compute_idu(a, x)
[m, n] = size(a); % size of truncation, a is mxn
ipi = intval('pi');
if ~isa(x,'taylor')
    k = size(x,2); % size of input, x is 1xk
    if m>1 % size of truncation greater than 1
        if k==1 % x is a number
            idbaseH1 = [1 sqrt(intval(2)).*cos((1:(m-1)).*ipi.*intval(x))];  % derivative of H^{1,0}_B base
            idu = idbaseH1*a; % 1xn vector function
        else % x is a vector
            A = intval(0).*zeros(k,m);% each row i of A is the base on x(i)
            for i=1:k
                A(i,:) = [1 sqrt(intval(2)).*cos((1:(m-1)).*ipi.*intval(x(i)))];
            end
            idu = reshape(A*a, [1,k*n]); % 1x(kxn), [A(x(1))*a(:,1),...,A(x(k))*a(:,1),...,A(x(1))*a(:,n),......,A(x(k))*a(:,n)]
        end
    else % size of truncation equal 1, that is, a is row vector
        idu = intval(0).*zeros(1,k*n); %1x(kxn) [a(1,1),...,a(1,1),...,(1,n),......,(1,n)]
        for i =1:n
            idu(((i-1)*k+1):i*k) = intval(a(1,i));
        end
    end
else
    v = x(1).t;
    taylorIndex = size(v,2)-1;
    f_aux = @(y) y;
    idu = f_aux(taylorinit(intval(1)*(1:n), taylorIndex));
    for j = 1:n
        idu(1,j) = intval(a(1,j));
        for i = 2:m
            idu(1,j) = idu(1,j) + sqrt(intval(2))*cos((i-1)*ipi*intval(x))*intval(a(i,j));
        end
    end
end
end
%%%%%%%%%%%%% vect_int_sympson function %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = vect_int_sympson(f, a, b, n)
m = size(f(a),2);
D = intval(b) - intval(a); H = D / n;
x = a + (intval(0:n)/n)*D;
w = 2 * ones(1,n+1); w(2:2:n) = 4; w(1) = 1; w(n+1) = 1;
vect_imgf = intval(zeros(n+1,m));
for i=1:(n+1)
    vect_imgf(i,:) = intval(f(x(i)));
end
V = H/3 .* (w * vect_imgf);
E = (H^4*D / intval('180'));
X = V- E; % Integral inclusion
end