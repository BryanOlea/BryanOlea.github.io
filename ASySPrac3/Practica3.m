%% Pr�ctica 3: Se�ales en tiempo discreto
% Caudillo Barbosa Eric
%
% Olea Garc�a Bryan
%
% Robles Mart�nez Dario Antonio
%% Objetivos
% 
% * Manipulaci�n b�sica de MATLAB
% * Gr�ficas de se�ales reales y complejas discretas
% * Transformaci�n de se�ales discretas (escalamientos y traslaciones)
% * Calculo de energ�a y potencia de se�ales discretas
% 
%% INTRODUCCI�N
% *Escalamiento horizontal de se�ales en tiempo discreto*
%%
% * *x2[n] = x1[an]*
%
% * Esta operaci�n originar� la aparici�n de nuevas muestras iguales a cero (a < 1) o la desaparici�n de algunas muestras (a > 1), debido a que la variable independiente "n" solo puede tomar valores enteros
% * Cuando se tiene la operaci�n de escalamiento en tiempo acompa�ada de un desplazamiento, primero se debe escalar la se�al y luego se debe desplazar. Estas operaciones tampoco son conmutativas entre s�.
%
% *Ejemplos:*
%
% *Si x(t) es una se�al de audio en una grabadora de cinta, x( 2 t ) ser�a la misma grabaci�n pero reproducida al doble de la velocidad) y x( � t ) reproducida a la mitad de la velocidad.
%
% *x2[n] = x1[ � n + 1]
n=-2:5;
x=[0 0 0 1 2 2 2 2 ];
Gd2(n,x)
title('x2[n]');
%%
n1=-6:8;
x1=[0 0 0 0 0 0 1 0 2 0 2 0 2 0 2];
Gd2(n1,x1)
title('x1[ � n + 1]');

%%
% *Desplazamiento*
%
% * *x2[n] = x1[n + n0]*
%
% * Equivale f�sicamente a adelantar o atrasar la se�al
% * Gr�ficamente equivale a desplazar la se�al hacia la izquierda (adelanto) o hacia la derecha (atraso).
% * El adelanto de una se�al no es posible f�sicamente, pero es muy �til su consideraci�n en el An�lisis de Se�ales
% * En la pr�ctica se pueden presentar dos casos: 
%
%  t0>0:Adelanto     t0<0:Atraso
%
% *Ejemplo:*
%
% x2[n] = x1[n + 2]
%
n2=-2:5;
x2=[0 0 0 1 2 2 2 2];
Gd2 (n2,x2)
title('x2[n]');
%%
n3=-2:5;
x3=[0 1 2 2 2 2 2 2];
Gd2 (n3,x3)
title('x1[n + 2]');
%%
% *Escalamiento en Magnitud*
% 
% *Equivale a multiplicar la se�al por una constante real.
% *En la pr�ctica se pueden presentar cuatro casos:
%
% A > 1 : Amplificador.
%
% A < 1 : Atenuador.
%
% A = 1 : Aislador.
%
% A = -1 : Inversor.
%
% *Ejemplo:*
%
% x2[n] = 2x1[n]
n4=-2:5;
x4=[0 0 0 1 2 2 2 2];
Gd2 (n4,x4)
title('x2[n]')
%%
n5=-2:5;
x5=[0 0 0 2 4 4 4 4];
Gd2 (n5,x5)
title('2x1[n]')
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