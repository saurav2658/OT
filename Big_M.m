clc
clear 
M=1000;
art_var=[5 6];
A=[1 3 -1 0 1 0; 1 1 0 -1 0 1];
b=[3; 2];
c=[-3 -5 0 0 -M -M 0]; 
a=[A b];
array2table(a,'VariableNames',{'x1','x2','s1','s2','A1','A2', 'b'});
bv=[5 6];
z=c(bv)*a-c;
simplex_table=[z;a];
Var={'x1','x2','s1','s2','A1','A2', 'b'};
array2table(simplex_table,'VariableNames',Var)

for k=1:15
if all(z(1:end-1)>=0)
   if any(bv==art_var(1))||any(bv==art_var(2))
       
       fprintf('Infeasible solution');
       break;
       
   end
    fprintf('The current table is optimal\n');
    optimal_value=z(end);
    fprintf('The optimal value of the current lpp is %f',optimal_value);
    break;
else
    fprintf('The current table is not optimal');
     [entering_var_value, pvt_col]=min(z(1:end-1));
     if all(a(:,pvt_col)<=0)
         error('The lpp is unbounded');
     else
         sol=a(:,end);
         column=a(:,pvt_col);
         for i=1:size(a,1)
             if column(i)>0
                 ratio(i)=sol(i)/column(i);
             else
                 ratio(i)=inf;
             end
         end
         [leaving_var_value,pvt_row]=min(ratio);
         bv(pvt_row)=pvt_col;
         pvt_key=a(pvt_row,pvt_col);
         a(pvt_row,:)=a(pvt_row,:)/pvt_key;
         for i=1:size(a,1)
             if i~=pvt_row 
                 a(i,:)=a(i,:)-a(i,pvt_col)*a(pvt_row,:);
             end
         end
         z=c(bv)*a-c;
         simplex_table=[z;a];
array2table(simplex_table,'VariableNames',Var)

     end
end
end

