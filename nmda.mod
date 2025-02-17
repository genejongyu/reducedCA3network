COMMENT
-----------------------------------------------------------------------------
Simple synaptic mechanism derived for first order kinetics of
binding of transmitter to postsynaptic receptors.

A. Destexhe & Z. Mainen, The Salk Institute, March 12, 1993.
-----------------------------------------------------------------------------

During the arrival of the presynaptic spike (detected by threshold 
crossing), it is assumed that there is a brief pulse (duration=Cdur)
of neurotransmitter C in the synaptic cleft (the maximal concentration
of C is Cmax).  Then, C is assumed to bind to a receptor Rc according 
to the following first-order kinetic scheme:

Rc + C ---(Alpha)--> Ro							(1)
       <--(Beta)--- 

where Rc and Ro are respectively the closed and open form of the 
postsynaptic receptor, Alpha and Beta are the forward and backward
rate constants.  If R represents the fraction of open gates Ro, 
then one can write the following kinetic equation:

dR/dt = Alpha * C * (1-R) - Beta * R					(2)

and the postsynaptic current is given by:

Isyn = gmax * R * (V-Erev)						(3)

where V is the postsynaptic potential, gmax is the maximal conductance 
of the synapse and Erev is the reversal potential.

If C is assumed to occur as a pulse in the synaptic cleft, such as

C     _____ . . . . . . Cmax
      |   |
 _____|   |______ . . . 0 
     t0   t1

then one can solve the kinetic equation exactly, instead of solving
one differential equation for the state variable and for each synapse, 
which would be greatly time consuming...  

Equation (2) can be solved as follows:

1. during the pulse (from t=t0 to t=t1), C = Cmax, which gives:

   R(t-t0) = Rinf + [ R(t0) - Rinf ] * exp (- (t-t0) / Rtau )		(4)

where 
   Rinf = Alpha * Cmax / (Alpha * Cmax + Beta) 
and
   Rtau = 1 / (Alpha * Cmax + Beta)

2. after the pulse (t>t1), C = 0, and one can write:

   R(t-t1) = R(t1) * exp (- Beta * (t-t1) )				(5)

There is a pointer called "pre" which must be set to the variable which
is supposed to trigger synaptic release.  This variable is usually the
presynaptic voltage but it can be the presynaptic calcium concentration, 
or other.  Prethresh is the value of the threshold at which the release is
initiated.

Once pre has crossed the threshold value given by Prethresh, a pulse
of C is generated for a duration of Cdur, and the synaptic conductances
are calculated accordingly to eqs (4-5).  Another event is not allowed to
occur for Deadtime milliseconds following after pre rises above threshold.

The user specifies the presynaptic location in hoc via the statement
	connect pre_GLU[i] , v.section(x)

where x is the arc length (0 - 1) along the presynaptic section (the currently
specified section), and i is the synapse number (Which is located at the
postsynaptic location in the usual way via
	postsynaptic_section {loc_GLU(i, x)}
Notice that loc_GLU() must be executed first since that function also
allocates space for the synapse.
-----------------------------------------------------------------------------
  GLUTAMATE SYNAPSE (AMPA-Kainate receptors)

  Parameters estimated from whole cell recordings of synaptic currents on
  Cochlear neurons (Raman & Trussel, Neuron 9: 173-186, 1992) as well as
  from sharp electrode EPSP's recordings in thalamocortical neurons (LGN)
  (Crunelli et al. J. Physiol. 384: 603, 1987).

-----------------------------------------------------------------------------
ENDCOMMENT



INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA
	POINTER pre
	RANGE C, R, R0, R1, g, gmax, lastrelease
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 1	(mM)		: max transmitter concentration
	Cdur	= 1.1	(ms)		: transmitter duration (rising phase)
	Alpha	= 10	(/ms mM)	: forward (binding) rate
	Beta	= 0.0125 (/ms)		: backward (unbinding) rate
	Erev	= 0	(mV)		: reversal potential
	Prethresh = 0 			: voltage level nec for release
	Deadtime = 0	(ms)		: mimimum time between release events
	gmax		(umho)		: maximum conductance
	eta     = 0.33  (/mM)
	mag     = 1     (mM)
	gamma   = 0.06  (/mV)
}

ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	R				: fraction of open channels
	R0				: open channels at start of release
	R1				: open channels at end of release
	Rinf				: steady state channels open
	Rtau		(ms)		: time constant of channel binding
	pre 				: pointer to presynaptic variable
	lastrelease	(ms)		: time of last spike
}

INITIAL {
	R = 0
	C = 0
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	lastrelease = -9e9
}

BREAKPOINT {
	SOLVE release
	g = (gmax * R)/(1 + eta * mag * exp( - (gamma * v)))
	i = g*(v - Erev)
}

PROCEDURE release() { LOCAL q
	:will crash if user hasn't set pre with the connect statement 

	q = ((t - lastrelease) - Cdur)		: time since last release ended

						: ready for another release?
	if (q > Deadtime) {
		if (pre > Prethresh) {		: spike occured?
			C = Cmax			: start new release
			R0 = R
			lastrelease = t
		}
						
	} else if (q < 0) {			: still releasing?
	
		: do nothing
	
	} else if (C == Cmax) {			: in dead time after release
		R1 = R
		C = 0.
	}



	if (C > 0) {				: transmitter being released?

	   R = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau)
				
	} else {				: no release occuring

  	   R = R1 * exptable (- Beta * (t - (lastrelease + Cdur)))
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}

FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000

	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}
    
    
