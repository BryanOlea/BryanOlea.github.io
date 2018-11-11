%% Pr�ctica 3: Se�ales en tiempo discreto
% Caudillo Barbosa Eric
% Olea Garc�a Bryan
%% Objetivos
% 
% * Manipulaci�n b�sica de MATLAB
% * Gr�ficas de se�ales reales y complejas discretas
% * Transformaci�n de se�ales discretas (escalamientos y traslaciones)
% * Calculo de energ�a y potencia de se�ales discretas
% 
%% Introducci�n
% DESCRIPTIVE TEXT
%% 1.-Crea una funci�n que se llame fun1 y reciba dos parametros $\Omega$ y $a$ la funci�n debe regresar la evaluaci�n $f[n]=r^ncos[\Omega n]+ r^nsen(\Omega n)j$, esta funci�n debe trabajar con $r\in R^+$, $\Omega\in R$ y $n\in  N^n$.
%   
%   function f= fun1( W,n,r )
%   f=(r.^n.*cos(W.*n))+(r.^n.*sin(W.*n).*i);
%   end
% 
%% 2.-Construya una funci�n que gr�fique funciones de $f:N\rightarrow R$ en el formato de su elecci�n y pruebe su c�digo mostrando la gr�fica de $x[n]=na^nu[n]$ vs $n$ para $n\in {-2,...,10}$ para $a=0.9$
f=@(n,a) (n.*a.^n).*heaviside(n);
n=-2:10;
Gd2(n,f(n,0.9))
%% 3.-Construya una funci�n que gr�fique funciones de $f:N\rightarrow R^2$ en el formato de su elecci�n y pruebe su c�digo mostrando la gr�fica de la funci�n exponencial del primer problema, no debe incluir el c�digo, solo el uso de la funci�n para mostrarla gr�fica. Reporte la gr�fica de $f[n]$ para $r=1.1$, $\Omega=0.5$ y $n\in{-2,...,20}$ (recuerde que ya tiene una funci�n para esto). Reporte la gr�fica de $|f[n]|$ vs $n$ y $\angle f[n]$ vs $n$
n=-2:20;
f=fun1(0.5,n,1.1);
Gd3(n,real(f),imag(f))
subplot(1,2,1)
stem(n,abs(f))
subplot(1,2,2)
stem(n,angle(f))
%% 4.-Programe una funci�n que calcule la energ�a de una se�al en tiemp discreto la fucion se debe llamar energiadis. La funci�n recibe dos param�tros de entrada: el vector de tiempo y las alturas asignadas. La funci�n regresa la energia de la se�al y despliega la gr�fica de la se�al.
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