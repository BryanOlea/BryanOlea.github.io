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
%% Crea una funci�n que se llame fun1 y reciba dos parametros $\Omega$ y $a$ la funci�n debe regresar la evaluaci�n $f[n]=r^ncos[\Omega n]+ r^nsen(\Omega n)j$, esta funci�n debe trabajar con $r\in R^+$, $\Omega\in R$ y $n\in  N^n$.
%   
%   function f= fun1( W,n,r )
%   f=(r.^n.*cos(W.*n))+(r.^n.*sin(W.*n).*i);
%   end
% 
%% Construya una funci�n que gr�fique funciones de $f:N\rightarrow R$ en el formato de su elecci�n y pruebe su c�digo mostrando la gr�fica de $x[n]=na^nu[n]$ vs $n$ para $n\in {-2,...,10}$ para $a=0.9$
f=@(n,a) (n.*a.^n).*heaviside(n);
n=-2:10;
Gd2(n,f(n,0.9))
%% Construya una funci�n que gr�fique funciones de $f:N\rightarrow R^2$ en el formato de su elecci�n y pruebe su c�digo mostrando la gr�fica de la funci�n exponencial del primer problema, no debe incluir el c�digo, solo el uso de la funci�n para mostrarla gr�fica. Reporte la gr�fica de $f[n]$ para $r=1.1$, $\Omega=0.5$ y $n\in{-2,...,20}$ (recuerde que ya tiene una funci�n para esto). Reporte la gr�fica de $|f[n]|$ vs $n$ y $\angle f[n]$ vs $n$
n=-2:20;
Gd3(n,real(fun1(0.5,n,1.1)),imag(fun1(0.5,n,1.1)))
subplot(1,2,1)
plot(n,abs(fun1(0.5,n,1.1)))
subplot(1,2,2)
plot(n,atan2(fun1(0.5,n,1.1)))
