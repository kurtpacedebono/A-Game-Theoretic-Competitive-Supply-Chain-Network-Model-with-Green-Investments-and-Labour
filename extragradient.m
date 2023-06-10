function [xstar xtmat]=extragradient(problem,zeta)
%Extragradient Algorithm for Solving VIs
eval(problem);

%Step 0: Initialisation
P=eye(n); %initialising identity matrix as input for quadprog
eps=1e-2; %setting tolerance
gap=0; %initialising gap
np=0; %initialising number of projections
options=optimset('LargeScale','off','Display','off'); %options for quadprog
xt_1 = x; %initial solution 
xtmat=[xt_1];%storing the iteration solutions in a matrix 
fxt_1=F; %initial function value 
nf=1; %initialising number of function evaluations 	
i=0; %initialising iteration counter 
tic %timing how long the algorithm takes
while  i==0 | norm(gap)>eps 
   i=i+1;
   %Step 1: Computation 
   F1=-(xt_1-zeta*fxt_1); %projection function
   [xbar]=quadprog(P,F1,A,b,Aeq,beq,lb,ub,xt_1,options); %solving projection
   xbar=abs(xbar);
   x=xbar;
   np=np+1; %iterating projections 
   if n==4
       fxbar=fs(x(1),x(2),x(3),x(4))';
   elseif n==11
       fxbar=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11))';
   elseif n==18
       fxbar=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18))';%evaluating function F(xbar) 
   elseif n==22
       fxbar=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18),x(19),x(20),x(21),x(22))';%evaluating function F(xbar) 
   end
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
   if n==4
       fxt_1=fs(x(1),x(2),x(3),x(4))';
   elseif n==11
       fxt_1=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11))';
   elseif n==18
       fxt_1=fs(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13),x(14),x(15),x(16),x(17),x(18))';%evaluating function F(xbar) 
   end
   nf=nf+1; %iterating function evaluation  
   fprintf('Iteration %g: |x^t - x^(t-1)|=%g, \n',i,norm(gap));%pause
   xtmat=[xtmat xt]; %storing solution
 end;
 
 %Tolerance reached: break out of loop
 fprintf('Number of Iterations=%g, Number of Projections=%g, Number of Function Evaluations=%g \n',i,np,nf)
 
 %Generate final solution 
 F=-(xt_1-fxt_1);
 [xstar]=quadprog(P,F,A,b,Aeq,beq,lb,ub,xt_1,options);
 xstar=abs(xstar);
 toc %returning the time taken
end