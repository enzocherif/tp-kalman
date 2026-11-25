

function y=k(x)

%.................DESCRIPTIF................
% Cette fonction est utilisée dans la plupart des blocs de Simulink
%
% Entrée: x: variable d'état qui représente la fraction molaire liquide
%
% Sortie: y:variable d'état qui représente la fraction molaire vapeur
%
%  La fonction k correspond à la loi non linéaire d'équilibre thermodynamique
%
%      sur [0 1], y=k(x)
%      sur ]-inf 0], on prolonge k(x) par y=k(0)+ k'(0)*x
%      sur [1 +inf[, on prolonge k(x) par y=k(1)+k'(1)*(x-1)


%.................CONSTANTE.................
% alpha: constante de volatilité relative
alpha=2.0;

%.................PROGRAMME.................
if and(x>=0,x<=1)
     y=(alpha*x)/(1+(alpha-1)*x);
   elseif x<0
     y=alpha*x;
   else
     y=1+(x-1)/alpha;
end;






