%% Práctica 3: Señales en tiempo discreto
% Caudillo Barbosa Eric
% Olea García Bryan
%% Objetivos
% 
% * Manipulación básica de MATLAB
% * Gráficas de señales reales y complejas discretas
% * Transformación de señales discretas (escalamientos y traslaciones)
% * Calculo de energía y potencia de señales discretas
% 
%% Introducción
% DESCRIPTIVE TEXT
%% 1.-Crea una función que se llame fun1 y reciba dos parametros $\Omega$ y $a$ la función debe regresar la evaluación $f[n]=r^ncos[\Omega n]+ r^nsen(\Omega n)j$, esta función debe trabajar con $r\in R^+$, $\Omega\in R$ y $n\in  N^n$.
%   
%   function f= fun1( W,n,r )
%   f=(r.^n.*cos(W.*n))+(r.^n.*sin(W.*n).*i);
%   end
% 
%% 2.-Construya una función que gráfique funciones de $f:N\rightarrow R$ en el formato de su elección y pruebe su código mostrando la gráfica de $x[n]=na^nu[n]$ vs $n$ para $n\in {-2,...,10}$ para $a=0.9$
f=@(n,a) (n.*a.^n).*heaviside(n);
n=-2:10;
Gd2(n,f(n,0.9))
%% 3.-Construya una función que gráfique funciones de $f:N\rightarrow R^2$ en el formato de su elección y pruebe su código mostrando la gráfica de la función exponencial del primer problema, no debe incluir el código, solo el uso de la función para mostrarla gráfica. Reporte la gráfica de $f[n]$ para $r=1.1$, $\Omega=0.5$ y $n\in{-2,...,20}$ (recuerde que ya tiene una función para esto). Reporte la gráfica de $|f[n]|$ vs $n$ y $\angle f[n]$ vs $n$
n=-2:20;
f=fun1(0.5,n,1.1);
Gd3(n,real(f),imag(f))
subplot(1,2,1)
stem(n,abs(f))
subplot(1,2,2)
stem(n,angle(f))
%% 4.-Programe una función que calcule la energía de una señal en tiemp discreto la fucion se debe llamar energiadis. La función recibe dos paramétros de entrada: el vector de tiempo y las alturas asignadas. La función regresa la energia de la señal y despliega la gráfica de la señal.
%
%   function E=energiadis(t,h)
%   E=sum(h.*h,2);
%   Gd2(t,h)
%   end
%
%% 5.-Resuelva el problema 3.1.1 c) de Lathi
t=-3:3;
n=[-9 -6 -3 0 3 6 9];
E=energiadis(t,n)
%% 6.-Resuelve el problema 3.2.3 de Lathi
x=@(n) n;
n=-10:10;
%-x[n]
Gd2(n,-x(n))
title('-x(n)','Interpreter','latex');
%x[n+6]
Gd2(n,x(n+6))
title('x(n+6)','Interpreter','latex');
%x[n-6]
Gd2(n,x(n-6))
title('x(n-6)','Interpreter','latex');
%x[3n]
Gd2(n,x(3.*n))
title('x(3n)','Interpreter','latex');
%x[n/3]
Gd2(n,x(n./3))
title('x(n/3)','Interpreter','latex');
%x[3-n]
Gd2(n,x(3-n))
title('x(3-n)','Interpreter','latex');