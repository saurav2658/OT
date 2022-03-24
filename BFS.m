clc
clear all
format short
A=[1 1 1 0;2 1 0 1]
c=[3 4 0 0]
b=[450;600]
n=size(A,2)
m=size(A,1)
if(n>m)
    ncm=nchoosek(n,m)
    pair=nchoosek(1:n,m)
    sol=[];
    for i=1:ncm
        y=zeros(n,1)
        x=A(:,pair(i,:))\b
        if all(x>=0 & x~=inf & x~=-inf)
            y(pair(i,:)) =x
            sol=[sol, y]
        end
    end
else
    error('ncm does not exists')
end
z=c*sol
[zmax, zindex]=max(z)
bfs=sol(:, zindex)
optimal_value=[bfs' zmax];
optimal_bfs=array2table(optimal_value)
optimal_bfs.Properties.VariableNames(1:size(optimal_bfs,2))={'x_1','x_2','x_3','x_4','z'}