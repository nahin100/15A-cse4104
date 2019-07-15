n1=[-1,0,1];
x1=[5,4,3];

n2=[-2,-1,0,1,2,3];
x2=[1,2,3,4,7,8,9];

m1 = min([n1,n2]);
m2 = max([n1,n2]);

m = m1:m2;

newx1= zeros(1, length(m));
newx2= zeros(1, length(m));
y = zeros(1, length(m));

for i = 1:length(m)
    for j = 1:length(n1)
        if m(i)==n1(j)
            newx1(i)=x1(j);
        end
    end
end

for i = 1:length(m)
    for j = 1:length(n2)
        if m(i)==n2(j)
            newx2(i)=x2(j);
        end
    end
end

y = newx1 + newx2;

n2=[-3,-2,-1,0,1,2,3];
x2=[1,2,3,4,7,8,9];

m = n2+3;

subplot(2,1,1); stem(n2,x2); xlabel('n'); ylabel('x(n)'); axis([-5,5,1,12]), set(gca, 'xtick', [-5 -4 -3 -2 -1 0 1 2 3 4 5]);
subplot(2,1,2); stem(m,x2); xlabel('n'); ylabel('x(n-3)'); axis([-5,5,1,12]), set(gca, 'xtick', [-5 -4 -3 -2 -1 0 1 2 3 4 5]);

%===============================

f = [1,2,3,1];
g = [1,2,1,-1];

% Transform the vectors f and g in new vectors with the same length
F = [f,zeros(1,length(g))];
G = [g,zeros(1,length(f))];
%C = zeros(1,length(F));

% FOR Loop to put the result of convolution between F and G vectors
% in a new vector C. According to the convolution operation characteristics,
% the length of a resultant vector of convolution operation between two vector
% is the sum of vectors length minus 1
for n=1:length(g)+length(f)-1
    % Create a new vector C
    C(n) = 0;
    % FOR Loop to walk through the vector F ang G
    for k=1:length(f)
        if(n-k+1>0)
            C(n) = C(n) + F(k) * G(n-k+1);
            fprintf("F(%d) = %d ", k, F(k));
            fprintf("G(%d) = %d\n",n-k+1, G(n-k+1));
            fprintf("C(%d) = %d\n",n, C(n));
        end
        fprintf("=======================\n");
    end
    fprintf("##########################\n");
    
end

% Show C vector on screen
C
 

%
        
        
        
    
