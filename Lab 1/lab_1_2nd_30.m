a=10;
b=20;
a+b
a/b
pi


a = [1, 2; 3, 4]
b = [1,2,3; 4,5,6]

row = size(b,1)
col = size(b,2)

b=2

row_vector = [1,2,3]
col_vector = [1;2;3]

b = [1,2,3; 4,5,6; 7,8,9]

% b(2,2) = 55
b(:,1)

%Scalar Matrix Operations

scalar = 2;
matrix = [1,4,7;3,6,9]

matrix-scalar

%Matrix Matrix Operations
m1 = [1,2;4,5]
m2 = [1,5,9;7,5,3]

mul = m1*m2

m1_transpose = m1'
m1

% Vector Creations

row_vector = [1,2,3]

zero_matrix = zeros(10,10)
one_matrix = ones(10,10)

rand_matrix = 10000*rand(10,10)

row_vector = 1:2:100

% Vector Operations
scalar = 2
row_vector = [1,2,3]

scalar*row_vector



rowv1 = [1,2,3];
rowv2 = [3,4,5];
rowv1.*rowv2

%Flip function & Rotate function

matrix = [1,2; 3,4; 5,6]

flip(matrix,2);
rot90(matrix)

%For Loop 

for i=1:2:10
    fprintf("i = %d\n", i);
end 

%While loop
    fprintf("========================\n");
i = 1;

sum=0;
mult=1;

% start = input('start = ');
% endd = input('end = ');
% increment = input('increment= ');
% 
% [sum,mult] = sum_and_mult(start, increment, endd);
%     
% fprintf("sum = %d\n", sum)
% fprintf("mult = %d\n", mult)



theta = 0:pi/100:2*pi;
y = sin(10*theta);
z = sin(10*theta)+cos(10*theta);

figure(1);
subplot(2,1,1), plot(x,y);
subplot(2,1,2), plot(x,z);


title('Two sine plot'),
xlabel('x data'),
ylabel('Amplitude');

set(gca, 'xlim', [0 2*pi]);
set(gca, 'ylim', [-3 +3]);

set(gca, 'xtick', [0.0 0.25 0.5 0.75 1.0]*2*pi);
set(gca, 'xticklabel', [{'0'} {'\pi/2'} {'\pi'} {'3\pi/2'} {'2\pi'}])

grid('on');
grid('minor');



% theta = 0:pi/100:2*pi;
% y = sin(20*theta);
% z = sin(20*theta)+cos(20*theta);
% 
% figure(2);
% plot(x, y, x,z);
% 
% title('Two sine plot'),
% xlabel('x data'),
% ylabel('Amplitude');
% 
% set(gca, 'xlim', [0 2*pi]);
% set(gca, 'ylim', [-3 +3]);
% 
% set(gca, 'xtick', [0.0 0.25 0.5 0.75 1.0]*2*pi);
% set(gca, 'xticklabel', [{'0'} {'\pi/2'} {'\pi'} {'3\pi/2'} {'2\pi'}])
% 
% grid('on');
% grid('minor');

a=10;
b=20000;
c=30;

if a>b && a>c
    fprintf("a = %d\n", a);
elseif b>a && b>c
    fprintf("b = %d\n", b);
else
    fprintf("c = %d\n", c);
end
































    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




















































































































