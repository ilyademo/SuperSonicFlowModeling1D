#ifndef GODUNOV_FUNCTIONS_H
#define GODUNOV_FUNCTIONS_H

extern "C" {
	void potok(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
	void raspad(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
	double divfgod(double*, double*, double*, double*);
	double fgod(double*, double*, double*, double*);
}

#endif // !GODUNOV_FUNCTIONS_H
