clc
clear
format rat

%Input parameter
%cost vector(c) 
c=[2,1];
%coeficient matrix(A) 
A=[1 2; 1 1; 1 -2];
B=[10;6;1];

%Plotting graph
p=max(B);
y1=0:1:max(B);
x11=(B(1)-A(1,1).*y1)./A(1,2);
x21=(B(2)-A(2,1).*y1)./A(2,2);
x31=(B(3)-A(3,1).*y1)./A(3,2);
x11=max(0,x11);
x21=max(0,x21);
x31=max(0,x31);
plot(y1,x11,'r',y1,x21,'b',y1,x31,'g')
title('x1 vs x2')
xlabel('value of x1')
ylabel('value of x2')

%Find corner point with axes
cx1=find(y1==0) %points with x1 axis
c1=find(x11==0) %points with x2 axis
Line1=[y1(:,[c1 cx1]); x11(:,[c1 cx1])]'

c2=find(x21==0) %points with x2 axis
Line2=[y1(:,[c2 cx1]); x21(:,[c2 cx1])]' ;

c3=find(x31==0) %points with x3 axis
Line3=[y1(:,[c3 cx1]); x31(:,[c3 cx1])]' ;
corpt=unique([Line1;Line2;Line3],'rows')

%Intersecting points
pt=[0;0];
for i=1:size(A,1)
    A1=A(i,:);
    B1=B(i,:);
for j=i+1:size(A,1)
     A2=A(j,:);
     B2=B(j,:);
     A4=[A1;A2]
     B4=[B1;B2]
     X=A4\B4
     pt=[pt X]
     end
     end
     ptt=pt'
     
 %All corner points
 allpt=[ptt;corpt]
 points=unique(allpt,'rows')
 
 %Feasible region
 PT=constraint(points)
 P=unique(PT,'rows')
 
 %Compute objective function
  %for i=1:size(P,1)
      %fn(i,:)=(sum(P(i,:).*C)%c multiplied by x then sum of them 
      %end
      %Vert_fns =[P fn];
 %Find the optimal one
 %[fxval ,indfx]= max(fn);
 %optval= Vert_fns(indfx,:);
 %Optimal_BFS= array2table(optval);
      