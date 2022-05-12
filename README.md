n_input = 2
m_input = 3
A_input = [3 1;4 3;1 2]
b_input = [3 6 4]
c_input = [-4 -1]
sign = [0 -1 1]
A = []
b = []
c = c_input
row_add = 0
for i = 1:m_input
     if sign(i) == 0
         A = [A;A_input(i,:)]
         b = [b b_input(i)]
         A = [A;-A_input(i,:)]
         b = [b -b_input(i)]
         row_add = row_add + 1
     else if sign(i) < 0
         A = [A;-A_input(i,:)]
         b = [b -b_input(i)]
     else
         A = [A;A_input(i,:)]
         b = [b b_input(i)]
     end
     end
end
m = m_input + row_add
n = n_input + m
S = eye(m)
A = [A S]
c = [c zeros(1,m)]
index = []
for i = 1:m
    index = [index n_input+i]
end
B = []
Cb = []
for i = 1:m
     B = [B A(:,index(i))];
     Cb = [Cb c(index(i))];
end
X = zeros(n,1);
Xb = inv(B)*b';
X(index) = Xb;
y = inv(B)*A;
Net_eval = Cb*y - c
for s = 1:100
     if Net_eval >= -10^(-10) 
         if Xb >= 0 
             disp("Solution is feasible")
             break
        else
             [x_min,LV] = min(Xb) 
             for j = 1:n
                 if y(LV,j) < 0
                    ratio(j) = abs(Net_eval(j)/y(LV,j))
                 else
                     ratio(j) = inf
                 end
            end
             [r_min,EV] = min(ratio) 
             if r_min == inf 
                 disp("Solution is infeasible")
                 break
             end
             index(LV) = EV
             B = []
             Cb = []
     for i = 1:m
         B = [B A(:,index(i))];
     Cb = [Cb c(index(i))];
     end
X = zeros(n,1);
Xb = inv(B)*b';
X(index) = Xb;
y = inv(B)*A;
Net_eval = Cb*y - c
        end
     else
         disp("Solution is not optimal")
         break
     end
end
disp(X)
c*X


Transportation Problem

clc
clear all
c = [2 10 4 5 ; 6 12 8 11 ; 3 9 5 7];
a = [ 12 25 20 ]
b = [ 25 10 15 5 ]
i = size(b,1);
j = size(a,1);
if sum(a) == sum(b)
    fprintf('Balanced');
else
    fprintf('Unbalanced');
        if(sum(b)<sum(a))
        c(:,end 1) = zeros(1,i)
        b(end 1) = sum(a) - sum(b)
    else
        c(end 1,:) = zeros(1,j)
        s(end 1) = sum(b) - sum(a)
        end
end
cc=c
m = size(c,1);
n = size(c,2);
x=zeros(m,n);
for i=1:m
    for j=1:n
cpq = min(min(c))
if cpq == inf
    break;
end
[p1,q1] = find(cpq==c)
xpq = min(a(p1),b(q1))
[x1,ind]=max(xpq)
p=p1(ind)
q=q1(ind)
x(p,q)=min(a(p),b(q))
if(x(p,q)==a(p))
    b(q)=b(q)-x(p,q)
    a(p)=a(p)-x(p,q)
    c(p,:)=Inf
else 
    b(q)=b(q)-x(p,q)
    a(p)=a(p)-x(p,q)
    c(:,q)=Inf
end
    end
end
c
x
z=x.*cc
cost=sum(sum(z))


Fibonacci

clc 
clear 
format short
f=@(x) (x<0.5).*((1-x)./2)+(x>=0.5).*(x.^2);
L=-1;
R=1;
n=6;
t=linspace(L,R,100);
plot(t,f(t));

fib=ones(1,n);
for i=3:n+1
    fib(i)=fib(i-1)+fib(i-2);
end
fib;

for k=1:n
    ratio=fib(n+1-k)./fib(n+2-k)
    x2=L+ratio*(R-L);
    x1=L+R-x2;
    fx1=f(x1)
    fx2=f(x2)
    rs1(k,:)=[ L R x1 x2 fx1 fx2]
    if fx1<fx2
        R=x2;
    elseif fx1>fx2
        L=x1;
    else if fx1==fx2
            if min(abs(x1),abs(L))==abs(L)
                R=x2;
            else
                L=x1;
            end
        end
    end
end

rsl(k,:)=[L R x1 x2 fx1 fx2]
var ={'L','R','x1','x2','fx1','fx2'};
res1=array2table(rs1);
res1.Properties.VariableNames(1:size(res1,2))=var;

xopt=(L+R)/2;
fopt=f(xopt);
xopt
fopt
