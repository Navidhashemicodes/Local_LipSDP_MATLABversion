clc
clear all
close all
yalmip('clear')



 F = [0.845182	0.238651;-0.477302	0.129228];
% F = [0.845182	0.238651;-2.1706 -0.5951];

pr=0.95;
r=sqrt(2*log(1/(1-pr)));
 L = [0.4812;-0.2936];
% L = [0.3899;-1.0092];


 Pp = [0.705984	-0.271997; -0.271997	0.320728];
% Pp = [2.1656 -2.5479; -2.5479 7.2887];




G = [0;1];
C = [1,0];

K = [1.6933 0.7243];

R1 = eye(2);
R2 = 1;

sgS = C*Pp*C' + R2;
sg = sqrt(sgS);

As=0.05;
sg = sqrt(sgS);
alpha = sqrt(2)*sg*erfinv(1-As);
kappa = alpha;



A=0.01:0.01:0.99;



for i=1:max(size(A))
    
clear p Pi sol SolverInfo LMI
a = A(i);
p = sdpvar(2,2);
% LMI=[a*p-F'*p*F -F'*p; -p*F -p+(((1-a)/(r^2))*eye(2))]>=0.000001*eye(4);
Navid=lmi([a*p-F'*p*F -F'*p; -p*F -p+(((1-a)/(r^2))*eye(2))]) 
% pe=[0.0164 0.0146;0.0146 0.0238];
% Navid=lmi([a*p-((F+G*K)')*p*(F+G*K)  ((F+G*K)')*p*G*K  -((F+G*K)')*p; K'*G'*p*(F+G*K)  ((1-a)*pe)-K'*G'*p*G*K  K'*G'*p; -p*(F+G*K)  p*G*K  -p]);
% Navid=addlmi(Navid,'[a*p-((F+G*K)'')*p*(F+G*K)  ((F+G*K)'')*p*G*K  -((F+G*K)'')*p; K''*G''*p*(F+G*K)            -K''*G''*p*G*K  K''*G''*p; -p*(F+G*K)  p*G*K  ((1-a)/(r^2))*eye(2)-p]>=0');
% LMI=[a*p-((F+G*K)')*p*(F+G*K)  ((F+G*K)')*p*G*K  -((F+G*K)')*p; K'*G'*p*(F+G*K)  ((1-a)*pe)-K'*G'*p*G*K  K'*G'*p; -p*(F+G*K)  p*G*K  -p]>=0.0000001*eye(6);
Navid=addlmi(Navid,'p>=0');
h=-log(det(p));
% sol=optimize(LMI,-logdet(p))
sol = solvesdp(Navid,h);
% sol = optimize(LMI,-geomean(p));


SolverInfo = sol.info;

P(:,:,i) = double(p);
Pi = inv(P(:,:,i));

if sol.problem == 1;
Vol(i) = 10000;    
else sol.problem == 0;
Vol(i) = det(Pi);
end


end


[Y,k]=min(Vol)
P(:,:,k)


% As = 0.10;

% Pe1 =    4.555614395277030   5.602951227252897
%          5.602951227252897   7.105896361583699

% Pe2 = 1.0e+02 *
% 
%    6.368263442855770   2.435201296689424
%    2.435201296689424   0.931779436203655





% As = 0.05;

% Pe1 =  3.822875057628987   4.701681378454622
%        4.701681378454622   5.962766598814131

%   Pe2 = 1.0e+02 *
% 
%    5.343593289788744   2.043368145167730
%    2.043368145167730   0.781851500803187






% As = 0.02;


% Pe1 =    3.220487713783610   3.960744550225317
%            3.960744550225317   5.022995943369243
% alpha = 3.0385


%    Pe2 = 1.0e+02 *
% 
%    4.501208768156635   1.721240891683949
%    1.721240891683949   0.658595278358821
%   alpha = 4.1391







