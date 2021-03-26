%%
%A Program to Calculate Impulse Response of
%a Blackhole on a Space-Body over the Time
%[Numerical Approach]
%%

%Put_Ins
NP = 100000;
tMax = 1;                               %[in s]
dt_Our= tMax/NP;                        %[in s]

c = 2.999*10^8;                         %[in m/s]
v0 = 10^6;                              %[in m/s]
G = 6.6743*10^-11;                      %[in m^3/(kg*s^2)]

m_S = 1.9889*10^30;                     %[in kg]
m_B = 100*m_S;                          %[in kg]
m_T = 555.55;                           %[in kg]

x_B = ((10^-15)/9.461)*10^1335;                             %[in m]
y_B = 0;                                %[in m]
z_B = 0;                                %[in m]

y_T = 0;                                %[in m]
z_T = 0;                                %[in m]
 

%Assigning_Arrays
t = zeros(1, NP);                       % A range from 0 to NP
Fg = zeros(1, NP);                      % A range from 0 to NP
x_T = zeros(1, NP);                     % A range from 0 to NP
r_BT = zeros(1, NP);                    % A range from 0 to NP
dt_T = zeros(1, NP);                    % A range from 0 to NP
v = zeros(1, NP);                       % A range from 0 to NP


%Initial_Conditions
x_T(1) = 0;                             %[in m]
v(1) = 0;                               %[in m/s]
t(1) = 0;                               %[in s]


%Looping
for i = 1: NP
    v(i + 1) = v(i) + 0.5*v0;                                                       %[in m/s]
  
    
    dt_T(i) = dt_Our*sqrt(1 - (v(i + 1)/c)^2)^0.5;                                  %[in s]
    
    
    x_T(i + 1) = x_T(i) + v(i + 1)*dt_T(i);                                         %[in m]
    
    
    r_BT(i + 1) = sqrt((x_B - x_T(i + 1))^2 + (y_B - y_T)^2 + (z_B - z_T)^2);       %[in m]
   
    
    Fg(i) = (G*m_B*m_T/(r_BT(i + 1))^2) + (m_T*(v(i + 1) - v(i))/dt_T(i));                                          
    
    
    t(i + 1) = t(i) + dt_T(i);                                                      %[in s]

end


%Plotting
plot(t(1: 100), Fg(1: 100), 'ro-');    %plotting X and Y
xlabel('Time [s]');                    %Labeling X
ylabel('Force_net [N]');               %Labeling Y
h = legend('Impulse Responce');        %Title of the plot
