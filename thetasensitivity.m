clear all
clc

syms x1 x2 x3 x4 v1 v2 v3 v4 v5 v6 v7 x5 x6 v8 v9 v10 v11 v12

%price fns
p11=100 - 0.001*(x1+x2)-0.0003*x5 +0.001*(v1+v2+v3+v4+v5+v6+v7)-0.0008*(v8+v9+v10+v11+v12);
p12=100 - 0.001*(x3+x4)-0.0003*x6 +0.0015*(v1+v2+v3+v4+v5+v6+v7)-0.001*(v8+v9+v10+v11+v12);

p21=130 - 0.0002*(x1+x2)-0.001*x5 -0.0008*(v1+v2+v3+v4+v5+v6+v7)+0.001*(v8+v9+v10+v11+v12);
p22=130 - 0.0002*(x3+x4)-0.001*x6 -0.001*(v1+v2+v3+v4+v5+v6+v7)+0.0015*(v8+v9+v10+v11+v12);

%cost fns
c1=0.00042*(x1+x3)^2+0.0005*(x1+x3)*v1;
c2=0.00032*(x2+x4)^2+0.0004*(x2+x4)*v2;
c3=0.00008*(x1+x3)^2+0.0001*(x1+x3)*v3;
c4=0.00008*(x2+x4)^2+0.0001*(x2+x4)*v4;
c5=0.0001*(x1+x2+x3+x4)^2+0.0005*(x1+x2+x3+x4)*v5;
c6=0.00008*(x1+x2)^2+0.0001*(x1+x2)*v6;
c7=0.00008*(x3+x4)^2+0.0001*(x3+x4)*v7;

c8=0.00052*(x5+x6)^2+0.0001*(x5+x6)*v8;
c9=0.0001*(x5+x6)^2+0.0001*(x5+x6)*v9;
c10=0.00009*(x5+x6)^2+0.0001*(x5+x6)*v10;
c11=0.00004*(x5)^2+0.0001*(x5)*v11;
c12=0.00004*(x6)^2+0.0001*(x6)*v12;

%labour cost fns
w1=(600/(75+0.01*v1))*(x1+x3);
w2=(600/(80+0.01*v2))*(x2+x4);
w3=(600/(300+0.01*v3))*(x1+x3);
w4=(600/(300+0.01*v4))*(x2+x4);
w5=(600/(150+0.01*v5))*(x1+x2+x3+x4);
w6=(600/(300+0.01*v6))*(x1+x2);
w7=(600/(300+0.01*v7))*(x3+x4);

w8=(700/(85+0.01*v8))*(x5+x6);
w9=(700/(450+0.01*v9))*(x5+x6);
w10=(700/(250+0.01*v10))*(x5+x6);
w11=(700/(450+0.01*v11))*(x5);
w12=(700/(450+0.01*v12))*(x6);

%emission tax fns
t1=5*(0.45*(x1+x3)-0.000005*v1*(x1+x3));
t2=5*(0.34*(x2+x4)-0.000004*v2*(x2+x4));
t3=5*(0.05*(x1+x3)-0.000001*v3*(x1+x3));
t4=5*(0.05*(x2+x4)-0.000001*v4*(x2+x4));
t5=5*(0.45*(x1+x2+x3+x4)-0.000005*v5*(x1+x2+x3+x4));
t6=5*(0.05*(x1+x2)-0.000001*v6*(x1+x2));
t7=5*(0.05*(x3+x4)-0.000001*v7*(x3+x4));

t8=5*(0.25*(x5+x6)-0.000003*v8*(x5+x6));
t9=5*(0.01*(x5+x6)-0.000001*v9*(x5+x6));
t10=5*(0.3*(x5+x6)-0.000004*v10*(x5+x6));
t11=0;
t12=0;

%profit fns
U1= p11*(x1+x2) + p12*(x3+x4) -c1-c2-c3-c4-c5-c6-c7-w1-w2-w3-w4-w5-w6-w7-t1-t2-t3-t4-t5-t5-t6-t7-(v1+v2+v3+v4+v5+v6);
U2= p21*(x5) + p22*(x6) -c8-c9-c10-c11-c12-w8-w9-w10-w11-w12-t8-t9-t10-t11-t12-(v8+v9+v10+v11+v12);

%defining F in VI(F,K) as the partial derivatives of U1 and U2 wrt the vars x and v
d_U1_x1=diff(-U1,x1);
d_U1_x2=diff(-U1,x2);
d_U1_x3=diff(-U1,x3);
d_U1_x4=diff(-U1,x4);
d_U1_v1=diff(-U1,v1);
d_U1_v2=diff(-U1,v2);
d_U1_v3=diff(-U1,v3);
d_U1_v4=diff(-U1,v4);
d_U1_v5=diff(-U1,v5);
d_U1_v6=diff(-U1,v6);
d_U1_v7=diff(-U1,v7);

d_U2_x1=diff(-U2,x5);
d_U2_x2=diff(-U2,x6);
d_U2_v8=diff(-U2,v8);
d_U2_v9=diff(-U2,v9);
d_U2_v10=diff(-U2,v10);
d_U2_v11=diff(-U2,v11);
d_U2_v12=diff(-U2,v12);


F_sym=[d_U1_x1, d_U1_x2, d_U1_x3,d_U1_x4,d_U1_v1,d_U1_v2,d_U1_v3,d_U1_v4,d_U1_v5,d_U1_v6,d_U1_v7,d_U2_x1,d_U2_x2,d_U2_v8,d_U2_v9,d_U2_v10,d_U2_v11,d_U2_v12];

theta=[];
thetasubs=0:5000:100000;

%dimension
n=18;
%note: matrices are in the form:
%[x1 x2 x3 x4, v1 v2 v3 v4 v5 v6 v7, x5 x6, v8 v9 v10 v11 v12]
%initial point
x0= [3000 3000 3000 3000, 0 0 0 0 0 0 0,3000 3000,0 0 0 0 0]';
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
    theta;
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
    10000;
    0;0;0;0;0];

b=[b1;b2];
Aeq = [];
beq = [];
%lower and upper bounds for x's and v's
lb = [0 0 0 0, 0 0 0 0 0 0 0,0 0,0 0 0 0 0];
ub = [];

xstarmat=[];

for j=1:size(thetasubs,2)
    theta=thetasubs(j);
        b = [9120; 9600; 36000; 36000; 18000; 36000;36000;
        theta;
        0;0;0;0;0;0;0;b2];
        %converting to anonymous fn for quadprog
        fs=matlabFunction(F_sym,'Vars',{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12});

        %Step 0: Initialisation
        P=eye(n); %initialising identity matrix as input for quadprog
        eps=1e-2; %setting tolerance
        gap=0; %initialising gap
        np=0; %initialising number of projections
        options=optimset('LargeScale','off','Display','off'); %options for quadprog 
        x=x0;
        xt_1=x0;
        fxt_1=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18))'; %initial function value 
        nf=1; %initialising number of function evaluations 	
        i=0; %initialising iteration counter 
        zeta=50; %inputting step size zeta
        tic %timing how long the algorithm takes
        while  i==0 | norm(gap)>eps 
           i=i+1;
           %Step 1: Computation 
           F1=-(xt_1-zeta*fxt_1); %projection function
           [xbar]=quadprog(P,F1,A,b,Aeq,beq,lb,ub,xt_1,options); %solving projection
           xbar=abs(xbar);
           x=xbar;
           np=np+1; %iterating projections 
           fxbar=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18))';%evaluating function F(xbar) 
           nf=nf+1; %iterating function evaluation 

           %Step 2: Adaptation
           F2=-xt_1+zeta*fxbar;%projection function
           [xt fval]=quadprog(P,F2,A,b,Aeq,beq,lb,ub,xt_1,options); %solving projection
           xt=abs(xt);
           np=np+1; %iterating projections 

           %Step 3: Convergence Verification 
           gap=xt-xt_1; %calculating gap
           xt_1=xt; %updating xt_1 with new value
           x=xt_1;
           fxt_1=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18))'; %updating fxt_1 with new value F(xt_1) 
           nf=nf+1; %iterating function evaluation  
         end;

         fprintf('Solution for theta= %g obtained\n',thetasubs(j));%pause
         F=-(xt_1-fxt_1);
         [xstar]=quadprog(P,F,A,b,Aeq,beq,lb,ub,xt_1,options);
         xstar=abs(xstar);
         xstarmat=[xstarmat'; xstar']';
end 
     
xstarmat

%%

fmat=[];
%Farm 1 emissions
e1=(0.45*(x1+x3)-0.000005*v1*(x1+x3))+(0.34*(x2+x4)-0.000004*v2*(x2+x4))+(0.05*(x1+x3)-0.000001*v3*(x1+x3))+(0.05*(x2+x4)-0.000001*v4*(x2+x4))+(0.45*(x1+x2+x3+x4)-0.000005*v5*(x1+x2+x3+x4))+(0.05*(x1+x2)-0.000001*v6*(x1+x2))+(0.05*(x3+x4)-0.000001*v7*(x3+x4));

%Farm 2 emissions
e2=(0.25*(x5+x6)-0.000003*v8*(x5+x6))+(0.01*(x5+x6)-0.000001*v9*(x5+x6))+(0.3*(x5+x6)-0.000004*v10*(x5+x6));

outmat=[];

for l=1:size(xstarmat,2)
    obj1=vpa(subs(U1,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(:,l)'))/100;
    obj2=vpa(subs(U2,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(:,l)'))/100;
    ps11=vpa(subs(p11,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(:,l)'))/100;
    ps12=vpa(subs(p12,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(:,l)'))/100;
    ps21=vpa(subs(p21,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(:,l)'))/100;
    ps22=vpa(subs(p22,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(:,l)'))/100;
    d11=xstarmat(1,l)+xstarmat(2,l);
    d12=xstarmat(3,l)+xstarmat(4,l);
    d21=xstarmat(12,l);
    d22=xstarmat(13,l);
    e1s=vpa(subs(e1,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(1:18,l)'));
    e2s=vpa(subs(e2,{x1,x2,x3,x4,v1,v2,v3,v4,v5,v6,v7,x5,x6,v8,v9,v10,v11,v12},xstarmat(1:18,l)'));
    output=[obj1 obj2 ps11 ps12 ps21 ps22 d11 d12 d21 d22 e1s e2s thetasubs(l)]';
    outmat=[outmat output];
    
    f1=xstarmat(1,l)+xstarmat(3,l);
    f2=xstarmat(2,l)+xstarmat(4,l);
    f3=xstarmat(1,l)+xstarmat(3,l);
    f4=xstarmat(2,l)+xstarmat(4,l);
    f5=xstarmat(1,l)+xstarmat(2,l)+xstarmat(3,l)+xstarmat(4,l);
    f6=xstarmat(1,l)+xstarmat(2,l);
    f7=xstarmat(3,l)+xstarmat(4,l);
    
    f8=xstarmat(12,l)+xstarmat(13,l);
    f9=
    f10=
    f11=
    f12=
    fmat=[fmat [f1 f2 f3 f4 f5 f6 f7]'];
end

roundn = @(x,n) round(x*10^n)./10^n;
outmatr=vpa(roundn(outmat,2))
%% 
xaxis=thetasubs;
plot(xaxis,max(outmatr(1,:)',0),'-o','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(2,:)','-s','LineWidth',1.5)
xlabel('Budget \theta')
ylabel('Profit')
legend('U1','U2')
%matlab2tikz('Theta1_profit.tex','showInfo', false);

%% 
figure 
plot(xaxis,outmatr(7,:)','-o','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(8,:)','-s','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(9,:)','-^','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(10,:)','-*','LineWidth',1.5)
xlabel('Budget \theta')
ylabel('Demand')
legend('D11','D12','D21','D22')
%matlab2tikz('Theta1_Demand.tex','showInfo', false);

%% 
figure 
plot(xaxis,outmatr(3,:)','-o','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(4,:)','-s','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(5,:)','-*','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(6,:)','-^','LineWidth',1.5)
xlabel('Budget \theta')
ylabel('Price')
legend('P11','P12','P21','P22')
%matlab2tikz('Theta1_Price.tex','showInfo', false);

%% 
figure 
plot(xaxis,outmatr(11,:)','-o','LineWidth',1.5)
hold on 
plot(xaxis,outmatr(12,:)','-s','LineWidth',1.5)
xlabel('Budget \theta')
ylabel('Total Emissions')
legend('Firm 1','Firm 2')
%matlab2tikz('Theta1_Emissions.tex','showInfo', false);

