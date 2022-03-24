clc
clear all
c=[7,4]
A=[1,-3;1,2;1,0]
B=[4;5;3]
S=eye(size(A,1))
I=[0,0,1]
index=find(I==1)
S(index,index)= -S(index,index)
mat=[A S B]
obj=array2table(c)
obj.Properties.VariableNames(1:size(c,2))={'x_1','x_2'}
cons=array2table(mat);
cons.Properties.VariableNames(1:size(mat,2))={'x_1','x_2','s1','s2','s3','B'}
