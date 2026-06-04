% Approaches to embryonic neurodevelopment: from neural cell to neural tube formation through mathematical models


% Fig.3. % How to convert lines to sphere manifolds  using bezier 

clc
clear
syms i j k v1 v2 v3 r s t rr ss tt W 

% The coordinates of the refrence vectors
 m1= [0,0,0];
 m2= [0,0,0];
 m3= [0,0,0];
 
 M1= [1,0,0];
 M2= [0,1,0];
 M3= [0,0,1];

 M11= [-1,0,0];
 M22= [0,-1,0];
 M33= [0,0,-1];
 rst1=[1];

 for rst=[-2.*pi:pi./24:2.*pi; 0:pi./48:2.*pi]

 % The primary sphere 
 X1P11= rst1.*cos(rst(1)).*sin(rst(2));
 Y1P21= rst1.*sin(rst(1)).*sin(rst(2));
 Z1P31= rst1.*cos(rst(2));
 %%%%%%%%%%%%
 X2P13= rst1.*cos(rst(1)+(2.*pi./3)).*sin(rst(2));
 Y2P23= rst1.*sin(rst(1)+(2.*pi./3)).*sin(rst(2));
 Z2P33= rst1.*cos(rst(2));
 %%%%%%%%%%%%%
 X3P13= rst1.*cos(rst(1)).*sin(rst(2));
 Y3P23= rst1.*sin(rst(1)).*sin(rst(2));
 Z3P33= rst1.*cos(rst(2));
 O=[0];

% Vector plot(M1,M2'-')
   quiver3(m1,m2,m3,M1,M2,M3,'AutoScale','off')
 hold on
   quiver3(m1,m2,m3,M11,M22,M33,'AutoScale','off')
 hold on


% Combinations: U=e1+e2+e3 & V=e11+e22+e33
% e1*e2 & e2*e3 & e1*e3 & e1*e33  
% Bezier point order:  e1*e0*e2 & e2*e0*e3 & e1*e0*e3, e1*e0*e33


tt= 0:.01:1;

% 1-  e1*e0* e2 
 %The W value could be also 1,-.5-.05i, 1, but we can change like
% 1,-.9-.05i, 1, for external buldging!
W(1)=[1];
W(2)=[-.3-.05i];
W(3)=[1];

hold on 
% 1-  e1*e0* e2


 R11=((X1P11).*W(1).*(tt - 1).^2 + (tt.^2).*(X2P13).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 R12=((Y1P21).*W(1).*(tt - 1).^2 + (tt.^2).*(Y2P23).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 R13=((Z1P31).*W(1).*(tt - 1).^2 + (tt.^2).*(Z2P33).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 plot3(R11,R12,R13)
hold on
  % 2-  e2*e0*e3

 R21=((X2P13).*W(1).*(tt - 1).^2 + (tt.^2).*(X3P13).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 R22=((Y2P23).*W(1).*(tt - 1).^2 + (tt.^2).*(Y3P23).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 R23=((Z2P33).*W(1).*(tt - 1).^2 + (tt.^2).*(Z3P33).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 plot3(R21,R22,R23)
hold on 
% 3-  e1*e0*e3

 R31=((X1P11).*W(1).*(tt - 1).^2 + (tt.^2).*(X3P13).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
 R32=((Y1P21).*W(1).*(tt - 1).^2 + (tt.^2).*(Y3P23).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
R33=((Z1P31).*W(1).*(tt - 1).^2 + (tt.^2).*(Z3P33).*W(3) - 2.*tt.*(O).*W(2).*(tt - 1))./((tt.^2).*W(3) + W(1).*(tt - 1).^2 - 2.*tt.*W(2).*(tt - 1));
  plot3(R31,R32,R33)
hold on 
%axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
 end

% Fig.4.B   %%%%%%%%%% Symmetry with zeta & Lie Group
clc
clear
x=-20:20
%x=sin(X)
%L=imag(exp(zeta(x)+x.*pi*.5*i))
for y=[-1:.1:1]
L=imag(zeta(exp((x).*pi.*y*i)))
hold on
plot(L)
xlabel('X')
ylabel('Y')
end

%%%%%%%%%%%% Fig.5.B
clc
clear
syms theta(t)
x = -.5:.1:.5;
t=[-i/2 -.4i -.3i -.2i -.1i 1i 2i 3i 4i i/2 .6i]
% C=-2.*pi:pi./2.7:2.*pi
% sig=0:pi./5:2.*pi
for y=[.1:.11:1] % y is 2 in the equation!!
T=(log(gamma((y*i*t+1)./4)-log(gamma((-y*i*t+1)./4)))./2i)-(log(pi*t)./2)  % Riemann_sigel Function  
% T=theta(t)
L=imag((exp((x).*pi.*y*i.*T)))
hold on
plot(x,L)
xlabel('X')
ylabel('Y')
end
L=imag((exp((x).*pi.*y*i.*T)))
hold on
plot(x,L)
xlabel('X')
ylabel('Y')

%%%%%%%%%%%%  Fig.5.C
clc
clear
x = -.5:.1:.5;
%n=[-i/2 -.4i -.3i -.2i -.1i 1i 2i 3i 4i i/2 .6i]
n=[-2 -4 -3 -2 -1 1 2 3 4 2 6]
t=(2.*pi.*(n+1-(7./8)))./(lambertw(1./exp(n+1-(7./8))))
% C=-2.*pi:pi./2.7:2.*pi
% sig=0:pi./5:2.*pi

for y=[.1:.11:1] % y is 2 in the equation!!      
T=(log(((y*i*t+1)./4)-log(((-y*i*t+1)./4)))./2i)-(log(pi*t)./2)  % Riemann_sigel Function  
L=complex((exp(i*T)).*zeta(.5+t*i))
hold on
plot(x,L)
xlabel('X')
ylabel('Y')
end

%%%%%%%%%%%% Figure.9
% Points on elliptic curve over finite field-2 with  neurons  

% Panel A
clc
clear
close all

disp('***** y^2 = x^3 + ax + b mod p *****') % p=15,a=3,b=0

p=input('insert prime p: ');
a=input('insert coefficient a: ');
b=input('insert coefficient b: ');

x=0:p-1;
y=0:p-1;

figure
hold on

for i=1:length(x)

    for j=1:length(y)

        if mod(y(j)^2 - x(i)^3 - a*x(i) - b,p)==0

            plot(x(i),y(j),...
                'ko',...
                'MarkerFaceColor','k',...
                'MarkerSize',8)

        end

    end

end

xlabel('x')
ylabel('y')

title(sprintf(...
'P(x,y)=y^2-x^3-%dx-%d \\equiv 0 (mod %d)',...
a,b,p))

axis equal
grid on


% The panel-B and B'
clc
clear

disp('***** Homogeneous Elliptic Curve *****');
disp('***** z*y^2 - x^3 - a*x*z^2 - b*z^3 = 0 (mod p) *****'); % p=16 a=3,b=2

n=input('insert prime p: ');
a=input('insert coefficient a: ');
b=input('insert coefficient b: ');

x=0:n-1;
y=0:n-1;
z=0:n-1;

figure
hold on

for i=1:n
    for j=1:n
        for k=1:n

            if mod( ...
                z(k)*(y(j)^2) ...
                - (x(i)^3) ...
                - a*x(i)*(z(k)^2) ...
                - b*(z(k)^3), ...
                n) == 0

                plot3(x(i),y(j),z(k),'.','MarkerSize',20)

            end
        end
    end
end

xlabel('x')
ylabel('y')
zlabel('z')

grid on
axis equal
view(3)


% Panel C
clc
clear
syms r s t
disp('***** Homogeneous Elliptic Curve *****');
disp('***** z*y^2 - x^3 - a*x*z^2 - b*z^3 = 0 (mod p) *****');

n=input('insert prime p: ');
a=input('insert coefficient a: ');
b=input('insert coefficient b: ');

x=[0:n-1];
y=[0:n-1];
z=[0:n-1];

figure
for i=1:n
    for j=1:n
        for k=1:n

            if mod( ...
                z(k)*(y(j)^2) ...
                - (x(i)^3) ...
                - a*x(i)*(z(k)^2) ...
                - b*(z(k)^3), ...
                n) == 0

                r=.7;

                [s,t]= meshgrid(-2.*pi:pi./20:2.*pi, -2.*pi:pi./20:2.*pi);

                ZX=((real((zeta(exp((r.*cos(s).*t.^5)./120 - (r.*cos(s).*t.^3)./6 + r.*cos(s).*t)).^2))).^(-1)).*(r.*cos(s).*sin(t));
                ZY=((real(((zeta(exp((r.*sin(s).*t.^5)./120 - (r.*sin(s).*t.^3)./6 + r.*sin(s).*t)).^2)))).^(-1)).*(r.*sin(s).*sin(t));
                ZZ=((real(((zeta(exp((r.*t.^4)./24 - (r.*t.^2)./2 + r)).^2)))).^(-1)).*(r.*cos(t));

                surf((ZX+x(i)),(ZY+y(j)),(ZZ+z(k)))
                hold on

            end
        end
    end
end

grid;
hold off;

% Panel C':Sphere / non-pyramidal neurons

clc
clear
syms r s t

disp('***** Homogeneous Elliptic Curve *****');
disp('***** z*y^2 - x^3 - a*x*z^2 - b*z^3 = 0 (mod p) *****');

n=input('insert prime p: ');
a=input('insert coefficient a: ');
b=input('insert coefficient b: ');

x=[0:n-1];
y=[0:n-1];
z=[0:n-1];

figure
for i=1:n
    for j=1:n
        for k=1:n

            if mod( ...
                z(k)*(y(j)^2) ...
                - (x(i)^3) ...
                - a*x(i)*(z(k)^2) ...
                - b*(z(k)^3), ...
                n) == 0

                r=1.9;

                [s,t]= meshgrid(-2.*pi:pi./20:2.*pi, -2.*pi:pi./20:2.*pi);

                ZX=-(sin(t)).*sin((real((zeta(exp((r.*cos(s).*t.^5)./120 - (r.*cos(s).*t.^3)./6 + r.*cos(s).*t)).^2))).^(-1)).*(r.*cos(s).*sin(t));
                ZY=-(sin(t)).*sin((real(((zeta(exp((r.*sin(s).*t.^5)./120 - (r.*sin(s).*t.^3)./6 + r.*sin(s).*t)).^2)))).^(-1)).*(r.*sin(s).*sin(t));
                ZZ=-(sin(t)).*sin((real(((zeta(exp((r.*t.^4)./24 - (r.*t.^2)./2 + r)).^2)))).^(-1)).*(r.*cos(t));

                surf((ZX+x(i)),(ZY+y(j)),(ZZ+z(k)))
                hold on

            end
        end
    end
end

grid;
hold off;
