%% Pr�ctica 5: Series de Fourier en tiempo continuo
% Caudillo Barbosa Eric
%
% Olea Garc�a Bryan
%
% Robles Mart�nez Dario Antonio
%% Objetivos
% * Realizar gr�ficas de series de Fourier exponenciales y trigonom�tricas en tiempo continuo
% * Manipulaci�n de instrucciones en MATLAB
% * Calculo n�merico de los coeficientes de Fourier
%% INTRODUCCI�N
%
% *Series de Fourier*
%
% �Que es la serie de Fourier?
%
% * Una funci�n $f(t)$ peri�dica de periodo $P$, se puede representar en forma de una suma infinita de funciones arm�nicas es decir
%
% $$f(t)=\frac{a_{0}}{b}+\sum_{k=1}^{\infty}(a_{k}cos(\frac{2\Pi
% }{P}t)+b_{k}sin(k\frac{2\Pi }{P}t))$$
%
% donde $a_{0} a_{1} ...a_{k} ...$ y $b_{1} b_{2} .... b_{k} ....$ son los denominados coeficientes de Fourier.
% (Una funci�n peri�dica, se puede representar en forma de una suma
% infinita de funciones arm�nicas).
%
% Si la funci�n f(t) tiene simetr�a, algunos de los coeficientes resultan nulos.
%
% * Si $f(t)$ es una funci�n par, $f(t)=f(-t)$, los t�rminos $b_{k}$ son nulos
%
% * Si $f(t)$ es impar $f(t)=-f(-t)$, los coeficientes $a_{k}$ son nulos
%%
%
% *Funci�n par*
%
% * Si la funci�n es par $b_{k}=0$
%
% Ejemplo:
%
% para el pulso rectangular sim�trico de anchura $1$ y periodo $P=2$ se obtienen los siguientes coeficientes.
%
figure
n=7; 
hold on
x=[-1 -0.5 -0.5 0.5 0.5 1];
y=[0 0 1 1 0 0];
plot(x,y,'b','linewidth',2)
x=linspace(-1,1,100);
y=zeros(length(x),1);
%%
%
% $$ a_{0}=\frac{2}{2}\int_{-.5}^{.5}dt=1 $$
%
% $$ a_{k}=\frac{2}{2}\int_{-.5}^{.5}cos(k\Pi t)dt=\frac{2}{k\pi}(sin(\frac{k\pi}{2}))=$$
%
% * $$ 0 k par$$
%
% * $$ \frac{2}{k\pi}(-1)^{(k-1)/2} k impar$$  
%
%%
syms t P k;
ak=int(cos(pi*k*t),t,-0.5,0.5);
subs(ak,k,sym('[1 2 3 4 5 6 7]'))
%%
%
% Vamos a reconstruir la funci�n f(t) a partir del desarrollo en serie de Fourier
%
% $$s_{n}(t)=\frac{1}{2}+\frac{2cos(t\pi)}{\pi}-\frac{2cos(3t\pi)}{3\pi}+\frac{2cos(5t\pi)}{5\pi}-\frac{2cos(7t\pi)}{7\pi}+...$$
%%
figure
n=7; 
hold on
x=[-1 -0.5 -0.5 0.5 0.5 1];
y=[0 0 1 1 0 0];
plot(x,y,'b','linewidth',2)
x=linspace(-1,1,100);
y=zeros(length(x),1);
for i=1:length(x)
    y(i)=1/2;
    for k=1:2:n
        y(i)=y(i)+(-1)^((k-1)/2)*2*cos(k*pi*x(i))/(k*pi);
     end
end
plot(x,y, 'r');
title(sprintf('Aproximaci�n de Fourier: %i t�rminos',n))
grid on
hold off
%%
% 
% *Funci�n impar*
%
% * Si la funci�n es impar,  $a_{k}=0$
%
% Sea ahora la funci�n de periodo $P=2$
%%
figure
n=7;
hold on
x=[-1 -1 0 0 1 1];
y=[0 1 1 -1 -1 0];
plot(x,y,'b','linewidth',2)
x=linspace(-1,1,100);
y=zeros(length(x),1);
%%
% Es una funci�n impar, los coeficientes a_{k} son nulos
% 
% $$b_{k}=\int_{-1}^{0}sin(k\pi t)dt-\int_{0}^{1}sin(k\pi
% t)dt=\frac{1}{k\pi}(-2+2cos(k\pi))=$$
%
% * $$ 0 k par$$
%
% * $$\frac{-4}{k\pi} k impar$$
%%
%
 syms t P k;
 bk=int(sin(pi*k*t),t,-1,0)-int(sin(pi*k*t),t,0,1);
 subs(bk,k,sym('[1 2 3 4 5 6 7]'))
%%
% El desarrollo en serie es:
%
% $$ s_{n}(t)=-\frac{4sin(t\pi)}{\pi}-\frac{4sin(3t\pi)}{3\pi}-\frac{4sin(5t\pi)}{5\pi}-\frac{4sin(7t\pi)}{7\pi}+... $$
% 
%%
%
figure
n=7;
hold on
x=[-1 -1 0 0 1 1];
y=[0 1 1 -1 -1 0];
plot(x,y,'b','linewidth',2)
x=linspace(-1,1,100);
y=zeros(length(x),1);
for i=1:length(x)
    y(i)=0;
    for k=1:2:n
        y(i)=y(i)-4*sin(k*pi*x(i))/(k*pi);
     end
end
plot(x,y, 'r');
title(sprintf('Aproximaci�n de Fourier: %i t�rminos',n))
grid on
hold off
%%
% *Aproximaci�n n�merica de los coeficientes de Fourier exponencial compleja*
%
% Podemos calcular $Dn$ num�ricamente mediante el uso de la DFT (la
% transformada discreta de Fourier)que utiliza las muestras de
% una se�al peri�dica $x (t)$ durante un per�odo. El intervalo de muestreo es de $T$ segundos.
% Por lo tanto, hay $$ N_{0}=\frac{T_{0}}{T} $$ n�mero de muestras en un per�odo
% $$ T_{0} $$. Para encontrar la relaci�n entre $$ D_{n} $$ y las muestras de $x
% (t)$
%
% $$ D_{n}=\frac{1}{T_{0}}\int_{T_{0}}^{0}x(t)e^{-jn\omega _{0}t}dt $$
%
% $$ =\lim_{T\rightarrow
% 0}\frac{1}{N_{0}T}\sum_{k=0}^{N_{0}-1}x(kT)e^{-jn\omega _{0}t}T $$
%
% $$ =\lim_{T\rightarrow
% 0}\frac{1}{N_{0}}\sum_{k=0}^{N_{0}-1}x(kT)e^{-jn\Omega _{0}k} $$
% 
% $$ \Omega _{0}=\omega _{0}T=\frac{2\pi}{N_{0}} $$
%
% En la pr�ctica, es imposible hacer que $$ T\rightarrow 0 $$ calcule el lado derecho de
% la ecuaci�n.Podemos hacer $T$ peque�a, pero no cero, lo que har�
% hacer que los datos aumenten sin l�mite. Por lo tanto, ignoraremos el
% l�mite de $T$ en la ecuaci�n  con el entendimiento impl�cito de que $T$ es
% razonablemente peque�o El valor distinto de cero $T$ dar� como resultado alg�n error de c�lculo, que es inevitable en cualquier evaluaci�n num�rica de una integral. los
% el error resultante de una $T$ distinta de cero se denomina error de alias
%
% $$ D_{n}=\frac{1}{N_{0}}\sum_{k=0}^{N_{0}-1}x(kT)e^{-jn\Omega _{0}k} $$
%%
T_0 = pi; N_0 = 256; T = T_0/N_0; t = (0:T:T*(N_0-1))'; M = 10;
x = exp(-t/2); x(1) = (exp(-pi/2) + 1)/2;
D_n = fft (x)/N_0; n = [-N_0/2:N_0/2-1]';
clf; subplot (2, 2, 1); stem(n, abs(fftshift (D_n)),'k');
axis ([-M M -.1 .6]); xlabel('n'); ylabel('|D_n|');
subplot (2, 2, 2); stem(n, angle(fftshift(D_n)),'k');
axis([-M M -pi pi]); xlabel ('n'); ylabel('\angle D n [rad]');
n = [0:M]; C_n(1) = abs(D_n(1)); C_n(2:M+1) = 2*abs (D_n(2:M+1));
theta_n(1) = angle(D_n(1)); theta_n(2:M+1) = angle(D_n(2:M+1));
subplot (2, 2, 3); stem(n,C_n,'k');
xlabel ('n'); ylabel('C_n');
subplot (2, 2, 4); stem(n,theta_n,'k');
xlabel ('n'); ylabel('\theta n [rad]');
%% 1.-Con serie y espectro trigonometrico, no es necesario entregar el c�digo, solo la aplicaci�n al problema especifico, debe de indicar la funci�n y los valores de sus coeficientes
% f(t)=e^{-t/2}; 0<t<\pi 

%% 2.-Con serie y espectro exponencial y A=3, no es necesario entregar el c�digo, solo la aplicaci�n al problema especifico, debe de indicar la funci�n y los valores de sus coeficientes
% f(t)=12triangularPulse(t-0.5)-6
% 4 arm�nicos
d0=-0.1
dn=@(n) ((-1.5/pi*n)*sin(0.5*pi*n)-(1.5/pi*n)*sin(1.5*pi*n)+(1.5*i/pi*n)*cos(0.5*pi*n)-(1.5*i/pi*n)*cos(1.5*pi*n)-(9*i/pi.^2*n.^2)*sin(0.5*pi*n)+(3*i/pi.^2*n.^2)*sin(1.5*pi*n)+(3/pi.^2*n.^2)*cos(0.5*pi*n)-(3/pi.^2*n.^2)*cos(1.5*pi*n));
t0=-0.5;
tf=1.5;
f=@(t) 12.*triangularPulse(t-0.5)-6;
armo=4;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)

% 15 arm�nicos
d0=-0.1
dn=@(n) ((-1.5/pi*n)*sin(0.5*pi*n)-(1.5/pi*n)*sin(1.5*pi*n)+(1.5*i/pi*n)*cos(0.5*pi*n)-(1.5*i/pi*n)*cos(1.5*pi*n)-(9*i/pi.^2*n.^2)*sin(0.5*pi*n)+(3*i/pi.^2*n.^2)*sin(1.5*pi*n)+(3/pi.^2*n.^2)*cos(0.5*pi*n)-(3/pi.^2*n.^2)*cos(1.5*pi*n));
t0=-0.5;
tf=1.5;
f=@(t) 12.*triangularPulse(t-0.5)-6;
armo=15;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)
%% 3.- Con serie y espectro exponencial, no es necesario entregar el c�digo, solo la aplicaci�n al problema especifico, debe de indicar la funci�n y los valores de sus coeficientes (sin incluir el procedimiento)
% f(t)=rectangularPulse(t/pi);
% 4 arm�nicos
d0=0.5
dn=@(n) (2*sin(pi*n/2)/(pi*n));
t0=-pi/2;
tf=pi/2;
f=@(t) rectangularPulse(t./pi);
armo=4;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)

% 15 arm�nicos
d0=0.5
dn=@(n) (2*sin(pi*n/2)/(pi*n));
t0=-pi/2;
tf=pi/2;
f=@(t) rectangularPulse(t./pi);
armo=15;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)
%% 4.-Con serie y espectro exponencial, no es necesario entregar el c�digo, solo la aplicaci�n al problema especifico, debe de indicar la funci�n y los valores de sus coeficientes (sin incluir el procedimiento)
% f(t)=|sin(t)|
% 4 arm�nicos
d0=2/pi;
dn=@(n) (2/pi*(1-4*n.^2));
t0=0;
tf=pi;
f=@(t) abs(sin(t));
armo=4;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)

% 15 arm�nicos
d0=2/pi;
dn=@(n) (2/pi*(1-4*n.^2));
t0=0;
tf=pi;
f=@(t) abs(sin(t));
armo=15;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)
%% 5.-Con serie y espectro exponencial y $T_0=3$ sin incluir la gr�fica de f, error ni energ�a del error, no es necesario entregar el c�digo, solo la aplicaci�n al problema especifico, debe de indicar la funci�n y los valores de sus coeficientes (sin incluir el procedimiento)
% f(t)=d(t)
% 4 arm�nicos
d0=1/3;
dn=1/3;
t0=0;
tf=3;
f=@(t) dirac(t);
armo=4;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)

% 15 arm�nicos
d0=1/3;
dn=1/3;
t0=0;
tf=3;
f=@(t) dirac(t);
armo=15;
a=-7;
b=7;
sfc(t0,tf,dn,d0,f,armo,a,b)
%% 6.-Elabore un c�digo similar al COMPUTER EXAMPLE C6.2 que se encuentra al final de la secci�n 6.2 de Lathi para el Ejempo 6.2 con los datos indicados anteriormente (no utilice inline)

x = [0 3 0 -3 0 3 0 -3 0];
tx = -2:0.5:2;
t = linspace (-2, 2,1000);
sumterms = zeros(16, length(t)); 
sumterms(1,:) = 0; %a0

for n = 1:size(sumterms,1)-1;
sumterms(n+1,:) = (((8*3)/((n^2)*(pi^2)))*sin(n*pi/2))*sin(n*pi*t);
end

x_N = cumsum (sumterms);
ind = 0;

for N = [0,1:2:size(sumterms, 1)-1],
ind = ind+1; 
subplot (3,3,ind);
plot (t,x_N(N+1,:),'b','linewidth',1);
hold on;
plot(tx,x,'k--');
axis ([-2 2 -4 4]);
grid on;
xlabel ('t'); 
ylabel (['x_{',num2str(N),'} (t)']);
end
%% 7.-Elabore un c�digo que implemente el algoritmo de trapecio compuesto para $n=15$, Utilice este c�digo para aproximar $D_0,...,D_4$ del ejemplo de la pr�ctica. Ahora implemente el c�digo COMPUTER EXAMPLE C6.4 que se encuentra al final de la secci�n 6.6 de Lathi, y calcule nuevamente el los coeficientes $D_0,...,D_4$ del ejemplo propuesto. Muestre una tabla que contenga los coeficientes mencionados calculados con los dos algoritmos y de forma exacta, �Qu� algortmo aproxima mejor a los coeficientes?, para esto compare los coefientes con el valor absoluto de la resta.
% clear all;
close all;
clc;

b=pi;
a=0;
n_lim = 15;
h = (b-a)/n_lim;
T0 = pi;

f = zeros(1,5);
w = 0:4;

%f(1,1)=0.504;

f1 = figure;

for n=0:4

f_0 = 1;
f_pi = exp(-2*j*n*pi)*exp(-pi/2);
fn = 0;


for t=1:n_lim-1
    fn = fn + exp(-2*j*n*(a+(h*t)))*exp(-((a+(h*t))/2));
end;

f(1,n+1) = (1/T0)*(h/2)*(f_0 + 2*fn + f_pi);

end;

w = 0:4;
subplot(2,1,1);
stem(w,real(f));
subplot(2,1,2);
stem(w,angle(f));


f2 = figure;

T_0 = pi; N_0 = 256; T = T_0/N_0; t = (0:T:T*(N_0-1))'; M = 10;
x = exp(-t/2); x(1) = (exp(-pi/2) + 1)/2;
D_n = fft (x)/N_0; n = [-N_0/2:N_0/2-1]';
subplot (2, 2, 1); stem(n, abs(fftshift (D_n)),'k');
axis ([-M M -.1 .6]); xlabel('n'); ylabel('|D_n|');
subplot (2, 2, 2); stem(n, angle(fftshift(D_n)),'k');
axis([-M M -pi pi]); xlabel ('n'); ylabel('\angle D n [rad]');
n = [0:M]; C_n(1) = abs(D_n(1)); C_n(2:M+1) = 2*abs (D_n(2:M+1));
theta_n(1) = angle(D_n(1)); theta_n(2:M+1) = angle(D_n(2:M+1));
subplot (2, 2, 3); stem(n,C_n,'k');
xlabel ('n'); ylabel('C_n');
subplot (2, 2, 4); stem(n,theta_n,'k');
xlabel ('n'); ylabel('\theta n [rad]');
