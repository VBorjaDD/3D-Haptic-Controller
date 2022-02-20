clear all;
close all;


addpath('/Users/louisjouret/Documents/GitHub/ADCS_LouisJouret/casadi-matlabR2015a-v3.3.0')
import casadi.*

model = stlread('Cellulo.stl');
Points = model.Points;
Points = Points/1000;
CList = model.ConnectivityList;


%% physical characteristics of the drone and constants
I1 = 0.00266;
I2 = 0.00266;
I3 = 0.00464;

kf = 1;
km = 1;
a = 0.2;
m = 1;
g = 9.81;

%% diff equations of the system

x = MX.sym('x',12);
para = MX.sym('para',9);

 ode =   [          x(4);
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

f = Function('f',{x,para},{ode},{'x','para'},{'ode'});


T   = 1.0; % Time horizon
N   = 10; % Number of control intervals
time= 10.0; %time simulated

% Integrator to discretize the system
intg_options = struct;
intg_options.tf = T/N; %integrate over T/N time
intg_options.number_of_finite_elements = 4; % Rugge-Kutta finite elements

% DAE problem structure
dae    = struct;
dae.x  = x;         % states
dae.p  = para;         % parameters (here the inputs)7
dae.ode= f(x,para);    % Expression for the right-hand side

intg   = integrator('intg','rk',dae,intg_options);
res    = intg('x0',x,'p',para); % Evaluate with symbols
x_next = res.xf;
% map the initial x&u to the next state
F      = Function('F',{x,para},{x_next},{'x','u'},{'x_next'});

V= VideoWriter('Drone2Movie','MPEG-4');
fps = N/T;
V.FrameRate = fps;

%% optimization loop

opti = casadi.Opti();
x    = opti.variable(12,N+1); % Decision variables for state trajectory
u    = opti.variable(4,N);
p    = opti.parameter(12,1);  % Parameter (not optimized over)
ref  = opti.parameter(12,1);
fm   = opti.parameter(5,N);

Q = diag([1 1 1 1 1 1 1 1 1 1 1 1]); 

para = [u; fm];
J = 0;

for k=1:N
    opti.subject_to(x(:,k+1)==F(x(:,k),para(:,k)));
    J = J + ((x(:,k)-ref).'*Q*(x(:,k)-ref));
end

opti.minimize(J); %minimize the cost function
opti.subject_to(x(:,1) == p); %first state is the previous state
opti.solver('ipopt');

%%
                
x_new =         [0; 0; 0; 0; 0; 0; 0; 0; 1.5; 0; 0; 0];
reference =     [0; 0; 0; 0; 0; 0; 0; 0; 1.5; 0; 0; 0];
Noise =         zeros(12,1);

X = x_new;
U = [];
Sol = [];
Ref = [];

ramp = [linspace(0,0,25*fps/10), linspace(0,1,50*fps/10), linspace(1,1,fps*time-50-25+N)];
force = [5*ramp; 0*ramp; -14*ramp];
M = [0.03*ramp; -0.03*ramp];
plot(linspace(0,fps*time,fps*time)/fps,force(:,1:fps*time)')

opti.set_value(p,x_new);
opti.set_value(ref, reference);
opti.set_value(fm,[force(:,1:N); M(:,1:N)]);

tspan = [0 1/fps];

for i=1:N*time/T
    sol     = opti.solve();
    input   = sol.value(u(:,1));
    Sol = [Sol sol.value(J)];
    inputs(1:4) = input;
    inputs(5:7) = force(:,i);
    inputs(8:9) = M(:,i);
    x_new   =  F(x_new,inputs);
    x_new   = full(x_new);
    X = [X x_new];
    U = [U input];
    opti.set_value(fm, [force(:,i:i+N-1);M(:,i:i+N-1)]);
    opti.set_value(ref, reference);
    opti.set_value(p,x_new); 
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
    axis(1*[-0.2,0.2,-0.2,0.2,1.2,1.6])
    frame = getframe;
    writeVideo(V,frame);
    i
end
close(V);