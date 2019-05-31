% IHC_demo.m
%
%						Inner Hair Cell Receptor Potential Model
%
%									David C. Mountain, Ph.D.
%										Boston University
%									Hearing Reseach Center
%
%  Model uses realistic transducer conductance and membrane time 
%  constant, but assumes basolateral membrane conductance is linear
%  This model assumes that the tension-gated channels are the only
%  apical channels.  Basolateral impedance is modeled as Rb/(Tau jw + 1)
%  and is implemented as a digital filter using the bilinear transform.
%
N = 2500;					% Number of time samples
Rate = 100000;				% Sampling rate (samples/s)
InArray  = zeros(1,N);	% Input waveform
Ga			= zeros(1,N);	% Apical conductance
OutArray = zeros(1,N);	% Output array for receptor potential
Time		= zeros(1,N);	% Time scale

% model parameters are from:
%
% Mountain, D.C and Cody, A.R. (1999)
% MULTIPLE MODES OF INNER HAIR CELL STIMULATION.
% Hearing Research 132: 1-14.

X0  = 27.0;		%{nm}
X1  = 27.0;		%{nm}
Sx0 = 85.0;		%{nm}
Sx1 = 11.0;		%{nm}
Gmax= 1.16e-8;		% 10 nS

Em	= -0.045;				% IHC basal membrane equilibrium potential = -45 mv
EP = .090-Em;				% driving potential = EP - IHC equilibrium potential = 125 mv
Rb = 1/58.8e-9;			% IHC basal resistance = 17 MOhms     
Tau = 0.0002;				% IHC time constant = 0.2 ms

Ga0   	= Gmax/(1+exp(X0/Sx0)*(1+exp(X1/Sx1)));	% apical conductance in quiet
deltaT 	= 1.0/Rate;			% Time step


A     	= 1/((Tau*2.0)/deltaT + 1);	%Coefficient used for membrane filter
B     	= Tau*2.0/deltaT - 1;			%Coefficient used for membrane filter

% Set up some initial conditions
Ga(1)   	= Ga0;
OldVm		= 0.0014;
Vm			= OldVm;
OldIR 	= 0.0014;
%
% Prompt the user for stimulus parameters
%
f = input('Frequency of tone (Hz) ==>');
amp=input('Cilia displacement (nm peak) ==>');
omega = 2*pi*f/1000;	%scale to kradians/s

Nstart 	= 500;	% start after delay in order to show resting potential
Nstop		= 2000;	% stop before end of record in order to show resting potential

Tstart	= deltaT*Nstart*1000;	% turn-on time for sinewave stimulus

for i = 1 : N					
   Time(i) = i*deltaT*1000;	% load the time array for plotting, scale to ms
	end;

for i = Nstart : Nstop						
   InArray(i) = amp*sin(omega*(Time(i)-Tstart));	% create input waveform
   end;
   
for i = 1 : N		%begin integration loop

   if i >1 OldIR = IR;
      end;
      Disp = InArray(i);
      Ga(i)   	= Gmax/(1+exp((X0-Disp)/Sx0)*(1+exp((X1-Disp)/Sx1)));
      IR = Ga(i)*(EP-Vm)*Rb;	%IR is the expected receptor potential if no membrane filtering
      OldVm = Vm;
      Vm    = A*(IR + OldIR+B*OldVm);			%this takes care of the membrane filtering
      OutArray(i) = (Em + Vm)*1000;				%scale from volts to mv and add Em

end;  %{integration loop}

Ga = 1e9*Ga;		% scale apical conductance to nS for plotting purposes

Fig=figure;				%Start generating the figure (assume 1024x768 display
set(Fig, 'Position',[5  100   600   600])	%Position it on the left side of screen
shg;

subplot(3,1,1)
plot(Time, InArray);
% autoscale axis since wide range of amplitudes can be used
grid;
ylabel('Displacment (nm)','FontSize',12,'FontWeight','bold');
title('IHC Model','FontSize',12,'FontWeight','bold');

subplot(3,1,2)
plot(Time, Ga);
axis([0 max(Time) 0 12]);
grid;
ylabel('Conductance (nS)','FontSize',12,'FontWeight','bold');

subplot(3,1,3)
plot(Time, OutArray);
axis([0 max(Time) -50 -20]);
grid;
xlabel('Time (seconds)','FontSize',12,'FontWeight','bold');
ylabel('Potential (mv)','FontSize',12,'FontWeight','bold');

