clear all;
close all;

time = 100;


V= VideoWriter('LQRMovie','MPEG-4');
fps = 10;
V.FrameRate = fps;


model = stlread('Cellulo concept 1 v47.stl');
Points = model.Points;
Points = Points/1000;
CList = model.ConnectivityList;

x = sym('x',[12 1]);
u = sym('u',[4 1]);

%% physical characteristics of the drone and constants
I1 = 4.5;
I2 = 4.3;
I3 = 6.0;

kf = 1;
km = 1;
a = 0.2;
m = 1;
g = 9.81;

%% small angle approx

ode =          [x(4);
                x(5);
                x(6);
                (1/I1)*(-x(6)*x(5)*(I2-I3) + kf*a*(u(2) - u(4)));
                (1/I2)*(-x(4)*x(6)*(I3-I1) + kf*a*(u(3) - u(1)));
                (1/I3)*(-x(4)*x(5)*(I1-I2) + km*(u(1) - u(2) + u(3) - u(4)));
                x(10);
                x(11);
                x(12);
                (kf/m)*(u(1) + u(2) + u(3) + u(4))*x(2);
                (kf/m)*(u(1) + u(2) + u(3) + u(4))*x(1);
                -g + (kf/m)*(u(1) + u(2) + u(3) + u(4))];
            
fun = @(t,x,u) [x(4);
                x(5);
                x(6);
                (1/I1)*(-x(6)*x(5)*(I2-I3) + kf*a*(u(2) - u(4)));
                (1/I2)*(-x(4)*x(6)*(I3-I1) + kf*a*(u(3) - u(1)));
                (1/I3)*(-x(4)*x(5)*(I1-I2) + km*(u(1) - u(2) + u(3) - u(4)));
                x(10);
                x(11);
                x(12);
                (kf/m)*(u(1) + u(2) + u(3) + u(4))*x(2);
                (kf/m)*(u(1) + u(2) + u(3) + u(4))*x(1);
                -g + (kf/m)*(u(1) + u(2) + u(3) + u(4))];


%% linearization for state-space representation
A = jacobian(ode, x)
B = jacobian(ode, u);
C = eye(12);
D = zeros(12,4);
h = 1/fps;

%% linearization around nominal value

A = jacobian(ode, x);
B = jacobian(ode, u);


h = 1/fps;

phi = h*A + eye(12);
gamma = h*B;

phi = subs(phi, x, [0;0;x(3);0;0;0;x(7); x(8); x(9);0;0;0]);
phi = double(subs(phi, u, [m*g/4;m*g/4;m*g/4;m*g/4]));
gamma = double(subs(gamma, x, [0;0;x(3);0;0;0;x(7); x(8); x(9);0;0;0]));


state =         [0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0];
reference =     [0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0];
NominalInput =  [m*g/4; m*g/4; m*g/4; m*g/4];
X = state;
U = NominalInput;
tspan = [0 h];


Q = diag([1,1,1,1,1,1,10,10,10,1,1,1]);

R = 1*diag([1,1,1,1]);

[K,S,e] = dlqr(phi,gamma,Q,R);

Ref = [reference];

for i = 1:fps*time
    input = -K*(state-reference) + NominalInput;
    input(input<0) = 0;
    input(input>20) = 20;
    [t,sol] = ode45(@(t, x) fun(t,x,input), tspan, state);
    oldref = reference;
    reference(7:9) = rotz(1)*reference(7:9);
    reference(9) = reference(9) + 0.01;
    reference(10:12) = (reference(7:9) - oldref(7:9))/h;
    Noise = 0.0001*(randn(12,1)-0.5);
    state = sol(end,:).' + Noise;
    X = [X state];
    U = [U input];
    Ref = [Ref reference];
    i
end



plot1 = figure;
plot(linspace(0,fps*time,fps*time + 1), X(1,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), X(2,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), X(3,:))
hold off;
xlabel('time [s]') 
ylabel('inclination [rad]') 
legend({'roll','pitch','yaw'})


plot2 = figure;
plot(linspace(0,fps*time,fps*time + 1), X(7,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), X(8,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), X(9,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1), Ref(7,:),'--')
hold on;
plot(linspace(0,fps*time,fps*time + 1), Ref(8,:),'--')
hold on;
plot(linspace(0,fps*time,fps*time + 1), Ref(9,:),'--')
hold off;
xlabel('time [s]') 
ylabel('coordinates [m]') 
legend({'x','y','z'})

% movie = figure;
% open(V);
% G = [];
% for i = 1:fps*time
%     Rx = rotx(X(1,i)*180/pi);
%     Ry = roty(X(2,i)*180/pi);
%     Rz = rotz(X(3,i)*180/pi);
%     rotatedPoints = Points*Rx*Ry*Rz;
%     G = [G  [X(7,i); X(8,i); X(9,i)]];
%     Translation = [X(7,i)*ones(length(Points),1), X(8,i)*ones(length(Points),1), X(9,i)*ones(length(Points),1)];
%     translatedPoints = rotatedPoints + Translation;
%     newPos = triangulation(CList,translatedPoints);
%     trisurf(newPos);
%     hold on;
%     fnplt(cscvn(G(:,1:end)),'r',2)
%     hold on;
%     scatter3(Ref(7,i),Ref(8,i),Ref(9,i),'*','g');
%     hold on;
%     fnplt(cscvn(Ref(7:9,1:i)),'g',2)
%     hold off;
%     axis(2*[-1,1,-1,1,-1,10])
%     frame = getframe;
%     writeVideo(V,frame);
%     i
% end
% close(V);