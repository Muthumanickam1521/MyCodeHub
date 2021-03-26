%-----------------------------------------------------------------------------------------------------------
% Numerical Computation on "Impulse Response of a Hypothetical Blackhole on a Body in Space-Time Linearly
%-----------------------------------------------------------------------------------------------------------



% Domain
NP = 100000;
tMax = 1;                              
dt_Our= tMax/NP;     


% Constants
c = 2.999*10^8;                         
v0 = 10^6;                              
G = 6.6743*10^-11;     

m_S = 1.9889*10^30;                     
m_B = 100*m_S;                          
m_b = 1000;         

x_B = ((10^-15)/9.461)*10^19;                             
y_B = 0;                               
z_B = 0;                                

y_b = 0;                               
z_b = 0;


% Assigning Arrays
t = zeros(1, NP);                    
Fg = zeros(1, NP);                      
x_b = zeros(1, NP);                    
r_Bb = zeros(1, NP);                    
dt_b = zeros(1, NP);                   
v = zeros(1, NP);                       


% Initial Conditions         
x_b(1) = 0;    
v(1) = 0;                               
t(1) = 0;                              


for i = 1: NP 
    v(i + 1) = v(i) + 0.5*v0;                                                     

    dt_b(i) = dt_Our*sqrt(1 - (v(i + 1)/c)^2);       

    x_b(i + 1) = x_b(i) + v(i + 1)*dt_b(i);                                        

    r_Bb(i + 1) = sqrt((x_B - x_b(i + 1))^2 + (y_B - y_b)^2 + (z_B - z_b)^2); 
    
    Fg(i) = (G*m_B*m_b/(r_Bb(i + 1))^2) + (m_b*(v(i + 1) - v(i))/dt_b(i));                                        

    t(i + 1) = t(i) + dt_b(i);        
end


%Plotting
plot(t(1: 100), Fg(1: 100), 'ro-');   
xlabel('Time [s]');                    
ylabel('Force net [N]');                
h = legend('Impulse Response');       
