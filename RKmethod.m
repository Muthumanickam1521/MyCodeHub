%------------------------------------------------
%Program to Apply 2-Step Range Kutta Method to
%Calculate Number of Nuclei in Radioactive Decay
%------------------------------------------------

%Put-Ins
N0 = 1000;
tau = 6;
NP = 100;
tMax = 10;

%Useful Relation
dt = tMax/NP;

%Assigning Arrays
N_Nuclei = zeros(1, NP);
t = zeros(1, NP);
t(1) = 0;

%Initialization
N_Nuclei(1) = N0;
t(1) = 0;

%Looping
for i = 1: NP
    m_1 = -N_Nuclei(i)/tau;
    m_2 = -(N_Nuclei(i) + (m_1)*(dt/2))/tau;
    N_Nuclei(i + 1) = N_Nuclei(i) + (m_2)*dt;
    t(i + 1) = t(i) + dt;
end

N_Analytics = N0*exp(-t/tau);

%Plotting
plot(t, N_Nuclei, 'go', t, N_Analytics, 'b');
xlabel('Time [s]');
ylabel('Nuclei Count');
h = legend('Time Vs Nuclei Count');
