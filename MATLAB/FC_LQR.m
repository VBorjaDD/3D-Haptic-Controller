clear all;
close all;

time = 15;


V= VideoWriter('LQRMovie','MPEG-4');
fps = 30;
V.FrameRate = fps;


model = stlread('Cellulo.stl');
Points = model.Points;
Points = Points/1000;
CList = model.ConnectivityList;

x = sym('x',[12 1]);
u = sym('u',[4 1]);
para = sym('para',[9,1]);

%% physical characteristics of the drone and constants
I1 = 0.00266;
I2 = 0.00266;
I3 = 0.00464;

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
                -(kf/m)*(u(1) + u(2) + u(3) + u(4))*x(1);
                -g + (kf/m)*(u(1) + u(2) + u(3) + u(4))];
            
fun = @(t,x,para) [ x(4);
                    x(5);
                    x(6);
                    (1/I1)*(-x(6)*x(5)*(I2-I3) + kf*a*(para(2) - para(4)) + para(8));
                    (1/I2)*(-x(4)*x(6)*(I3-I1) + kf*a*(para(3) - para(1)) + para(9));
                    (1/I3)*(-x(4)*x(5)*(I1-I2) + km*(para(1) - para(2) + para(3) - para(4)));
                    x(10);
                    x(11);
                    x(12);
                    (kf/m)*(para(1) + para(2) + para(3) + para(4))*sin(x(2)) + para(5)/m;
                    -(kf/m)*(para(1) + para(2) + para(3) + para(4))*cos(x(2))*sin(x(1)) + para(6)/m;
                    -g + (kf/m)*(para(1) + para(2) + para(3) + para(4))*cos(x(1))*cos(x(2)) + para(7)/m];


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


state =   [0; 0; 0; 0; 0; 0; 0; 0; 1.5; 0; 0; 0];
reference = state;
NominalInput =  [m*g/4; m*g/4; m*g/4; m*g/4];
X = state;
U = [];
tspan = [0 h];


Q = diag(1*ones(1,12));

R = 0*diag([1,1,1,1]);

[K,S,e] = dlqr(phi,gamma,Q,R);
intg = zeros(12,1);

inputs = zeros(9,1);

ramp = [linspace(0,0,25), linspace(0,1,50), linspace(1,1,time*fps-50-25)];
force = [5*ramp; 0*ramp; -14*ramp];
plot(linspace(0,fps*time,fps*time)/fps,force')
xlabel('time [s]') 
ylabel('Force [N]') 
legend({'Fx','Fy','Fz'})

M = [0.03*ramp; -0.03*ramp];


for i = 1:fps*time
    inputIntg = - 0.05*K*intg;
    input = -K*(state -reference) + NominalInput + inputIntg;
    input(input<0) = 0;
    inputs(1:4) = input;
    inputs(5:7) = force(:,i);
    inputs(8:9) = M(:,i);
    [t,sol] = ode45(@(t, x) fun(t,x,inputs), tspan, state);
    Noise = 0.00*(randn(12,1)-0.5);
    intg(7:9) = intg(7:9) + (state(7:9) -reference(7:9));
    state = sol(end,:).' + Noise;
    X = [X state];
    U = [U input];
    i
end



plot1 = figure;
plot(linspace(0,fps*time,fps*time + 1)/fps, X(1,:)*180/pi)
hold on;
plot(linspace(0,fps*time,fps*time + 1)/fps, X(2,:)*180/pi)
hold on;
plot(linspace(0,fps*time,fps*time + 1)/fps, X(3,:)*180/pi)
hold off;
xlabel('time [s]') 
ylabel('inclination [deg]') 
legend({'roll','pitch','yaw'})

plot4 = figure;
plot(linspace(0,fps*time,fps*time)/fps, U(1,:))
hold on;
plot(linspace(0,fps*time,fps*time)/fps, U(2,:))
hold on;
plot(linspace(0,fps*time,fps*time)/fps, U(3,:))
hold on;
plot(linspace(0,fps*time,fps*time)/fps, U(4,:))
hold off;
xlabel('time [s]') 
ylabel('input') 
legend({'1','2','3','4'})

plot2 = figure;
plot(linspace(0,fps*time,fps*time + 1)/fps, X(7,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1)/fps, X(8,:))
hold on;
plot(linspace(0,fps*time,fps*time + 1)/fps, X(9,:))
hold off;
xlabel('time [s]') 
ylabel('coordinates [m]') 
legend({'x','y','z'})


movie = figure;
open(V);
G = [];
for i = 1:fps*time
    Rx = rotx(X(1,i)*180/pi);
    Ry = roty(X(2,i)*180/pi);
    Rz = rotz(X(3,i)*180/pi);
    rotatedPoints = Points*Rx*Ry*Rz;
    G = [G  [X(7,i); X(8,i); X(9,i)]];
    Translation = [X(7,i)*ones(length(Points),1), X(8,i)*ones(length(Points),1), X(9,i)*ones(length(Points),1)];
    translatedPoints = rotatedPoints + Translation;
    newPos = triangulation(CList,translatedPoints);
    trisurf(newPos);
    hold on;
    fnplt(cscvn(G(:,1:end)),'r',2)
    hold off;
    axis(1*[-1,1,-1,1,0,2])
    frame = getframe;
    writeVideo(V,frame);
    i
end
close(V);
