clc
clear all

a=[1 1;3 2]
B=[5;12]
c=[6 5]
Noofvariables=2
s=eye(size(a,1))
A=[a s B]
cost=zeros(1,size(A,2))
cost(1:Noofvariables)=c
bv=Noofvariables+1:1:size(A,2)-1
zjcj=cost(bv)*A-cost
zcj=[zjcj;A];
simptable=array2table(zcj)
simptable.Properties.VariableNames(1:size(zcj,2))={'x_1' ,'x_2' ,'s_1','s_2','sol'}
RUN=true
while RUN
zc=zjcj(1:end-1);
if any(zc<0);
    fprintf('the current BFS is not optimal \n')
    
    [Enter_val, pvt_col]= min(zc)
    if all(A(:,pvt_col)<0)
        error('LPP is unbounded all entries are <=0 in column %d',pvt_col);
    else
        sol=A(:,end)
        column=A(:,pvt_col)
        for i=1:size(A,1)
            if column(i)>0
                ratio(i)= sol(i)./column(i)
            else
                ratio(i)=inf
            end
        end
        [leaving_val,pvt_row]=min(ratio)
    end
    bv(pvt_row)=pvt_col;
    pvt_key=A(pvt_row, pvt_col);
    A(pvt_row,:)=A(pvt_row,:)./pvt_key
    for i=1:size(A,1)
        if i~=pvt_row
            A(i,:)=A(i,:)-A(i, pvt_col).*A(pvt_row,:);
        end
    end
    zjcj=zjcj-zjcj(pvt_col).*A(pvt_row,:);
    zcj=[zjcj;A];
    table=array2table(zcj);
    table.Properties.VariableNames(1:size(zcj,2))={'x_1','x_2','s_1','s_2','sol'}
else
    RUN=false;
    fprintf('The current BFS is optimal \n')
end
end