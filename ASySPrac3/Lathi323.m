%% 6.-Resuelve el problema 3.2.3 de Lathi
% *a) x[-n]*
%%
% *b) x[n+6]*
%%
% *c) x[n-6]*
%%
% *d) x[3n]*
%%
% *e) x[n/3]*

n=0:6;
x=[0 1 2 3 2 1 0];
Gd2(n,x)
title('x[n]');
%%
Gd2(-n,x)
title('x[-n]');
%%
Gd2(n-6,x)
title('x[n+6]');
%%
Gd2(n+6,x)
title('x[n-6]');
%%
x1=[0 3 0 0 0 0 0];
Gd2(n,x1)
title('x[3n]');
%%
Gd2(3*n,x)
title('x[n/3]');