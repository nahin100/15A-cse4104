function [sum, mult] = sum_and_mult(start,increment,endd)

if (nargin~=3)
    disp('error');
end

sum=0;
mult=1;

for i = start:increment:endd
    sum=sum+i;
    mult=mult*i;
end

