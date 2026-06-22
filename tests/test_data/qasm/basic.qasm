OPENQASM 3.0;
#PRAGMA I <3 QGL

qubit q1;
length t;

//A very special gate
gate X(angle[32]: phi) q {
	U(phi, -pi/2, -pi/2) q;
}

/* Let's do a ramsey experiment!
   Fun! */

for t in 4ns:10us:20ns {
	reset    q1;
	X(pi/2)  q1;
	delay[t] q1;
	X(pi/2)  q1;
	measure  q1;
}