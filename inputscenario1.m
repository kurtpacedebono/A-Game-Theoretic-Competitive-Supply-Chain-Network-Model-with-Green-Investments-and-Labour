%dimension
n=18;
%note: matrices are in the form:
%[x1 x2 x3 x4, v1 v2 v3 v4 v5 v6 v7, x5 x6, v8 v9 v10 v11 v12]
%initial point
x= [3000 3000 3000 3000, 0 0 0 0 0 0 0,3000 3000,0 0 0 0 0]';
%FARM 1: 
%flow constraints 
A(1,:)=[1 0 1 0, -1.2 0 0 0 0 0 0,0 0, 0 0 0 0 0];
A(2,:)=[0 1 0 1, 0 -1.2 0 0 0 0 0,0 0, 0 0 0 0 0];
A(3,:)=[1 0 1 0, 0 0 -1.2 0 0 0 0,0 0, 0 0 0 0 0];
A(4,:)=[0 1 0 1, 0 0 0 -1.2 0 0 0,0 0, 0 0 0 0 0];
A(5,:)=[1 1 1 1, 0 0 0 0 -1.2 0 0,0 0, 0 0 0 0 0];
A(6,:)=[1 1 0 0, 0 0 0 0 0 -1.2 0,0 0, 0 0 0 0 0];
A(7,:)=[0 0 1 1, 0 0 0 0 0 0 -1.2,0 0, 0 0 0 0 0];

%budget inequality constraint 
A(8,:)=[0 0 0 0, 1 1 1 1 1 1 1,0 0, 0 0 0 0 0];

%vmax constraints  
A(9,:)=[-0.5 0 -0.5 0  ,1 0 0 0 0 0 0,0 0, 0 0 0 0 0];
A(10,:)=[0 -0.5 0 -0.5  ,0 1 0 0 0 0 0,0 0, 0 0 0 0 0];
A(11,:)=[-0.5 0 -0.5 0  ,0 0 1 0 0 0 0,0 0, 0 0 0 0 0];
A(12,:)=[0 -0.5 0 -0.5  ,0 0 0 1 0 0 0,0 0, 0 0 0 0 0];
A(13,:)=[-0.5 -0.5 -0.5 -0.5,0 0 0 0 1 0 0,0 0, 0 0 0 0 0];
A(14,:)=[-0.5 -0.5 0 0  ,0 0 0 0 0 1 0,0 0, 0 0 0 0 0];
A(15,:)=[0 0 -0.5 -0.5  ,0 0 0 0 0 0 1,0 0, 0 0 0 0 0]; 


b1 = [9120; 9600; 36000; 36000; 18000; 36000;36000;
    5000;
    0;0;0;0;0;0;0];


%FARM 2: 
%flow constraints 
A(16,:)=[0 0 0 0, 0 0 0 0 0 0 0,1 1,-2.5 0 0 0 0];
A(17,:)=[0 0 0 0, 0 0 0 0 0 0 0,1 1,0 -2.5 0 0 0];
A(18,:)=[0 0 0 0, 0 0 0 0 0 0 0,1 1,0 0 -2.5 0 0];
A(19,:)=[0 0 0 0, 0 0 0 0 0 0 0,1 0,0 0 0 -2.5 0];
A(20,:)=[0 0 0 0, 0 0 0 0 0 0 0,0 1,0 0 0 0 -2.5];
    
%budget inequality constraint 
A(21,:)=[0 0 0 0, 0 0 0 0 0 0 0,0 0,1 1 1 1 1];

%vmax constraints  
A(22,:)=[0 0 0 0, 0 0 0 0 0 0 0,-0.5 -0.5,1 0 0 0 0];
A(23,:)=[0 0 0 0, 0 0 0 0 0 0 0,-0.5 -0.5,0 1 0 0 0];
A(24,:)=[0 0 0 0, 0 0 0 0 0 0 0,-0.5 -0.5,0 0 1 0 0];
A(25,:)=[0 0 0 0, 0 0 0 0 0 0 0,0 0,0 0 0 1 0];
A(26,:)=[0 0 0 0, 0 0 0 0 0 0 0,0 0,0 0 0 0 1];

b2 = [21250; 112500; 62500;112500;112500;
    5000;
    0;0;0;0;0];

b=[b1;b2];
Aeq = [];
beq = [];
%lower and upper bounds for x's and v's
lb = [0 0 0 0, 0 0 0 0 0 0 0,0 0,0 0 0 0 0];
ub = [];
%function F from VI(F,K) 
fs = eval('scenario1');
F = fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18))';