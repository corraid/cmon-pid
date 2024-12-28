#pragma once

// Copyright 2024 corraid
//
// Redistribution and use in source and binary forms, with or without modification, are 
// permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of 
// conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of 
// conditions and the following disclaimer in the documentation and/or other materials 
// provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// PID transfer function realization with backward euler integration.
class pid_bwe {
protected:
	double D;        // lowpass filter
	double I;        // integrator
	double A1, A3, B3, C3;
public:
	// Update the pid with a new input value 'e'.
	// It advances the pid state by the time step 'h'.
	// When used in real time this function must be called every time 'h' has elapsed.
	// The return value is the control output 'u' of the transfer function u(s)/e(s).
	double Update(double e) {
		D = A3 * e + A1 * D;
		I = B3 * e + I;
		return C3 * e + I + D;
	};
	
	// Setup parameters for a serial pid with the transfer function
	// u(s) / e(s) = G * ((Ti * s + 1) / (Ti * s) * (Td * s + 1) * (1 / (Tf * s + 1)));
	// h is the sampling period. It has the same unit as Ti, Td and Tf.
	// G is the pid gain.
	void SerialPid(double h, double G, double Ti, double Td, double Tf) {
		A1 = Tf / (h + Tf);
		A3 = (G * (Tf - Td) * (Ti - Tf) * h) / (Tf * Ti * (h + Tf));
		B3 = (G * h) / Ti;
		C3 = (G * Td) / Tf;
	}
	
	// Setup parameters for a parallel pid with the transfer function
	// u(s) / e(s) = (Kp + Ki / s + Kd * s) * 1 / (Tf * s + 1)
	// h is the sampling period.
	// Tf is the filter time constant.
	void ParallelPid(double h, double Kp, double Ki, double Kd, double Tf) {
		A1 = Tf / (h + Tf);
		A3 = -(((Ki * Tf - Kp) * Tf + Kd) * h) / (Tf * (h + Tf));
		B3 = Ki * h;
		C3 = Kd / Tf;
	}

	// Setup parameters for a standard pid with the transfer function
	// u(s) / e(s) = Kp * (1 + 1 / (Ti * s) + s * Td / ((Td / N) * s + 1))
	// h is the sampling period. It has the same unit as Ti and Td.
	// Td/N is the filter time constant.
	// Kp is the pid gain.
	void NStandardPid(double h, double Kp, double Ti, double Td, double N) {
		A1 = Td / (N * h + Td);
		A3 = -(Kp * N * N * h) / (N * h + Td);
		B3 = (Kp * h) / Ti;
		C3 = Kp * N + Kp;
	}

	// Setup parameters for a standard pid with the transfer function
	// u(s) / e(s) = Kp * (1 + 1 / (Ti * s) + s * Td) * (1 / (Tf * s + 1))
	// h is the sampling period. It has the same unit as Ti, Td and Tf.
	// Tf is the filter time constant.
	// Kp is the pid gain.
	void StandardPid(double h, double Kp, double Ti, double Td, double Tf) {
		A1 = Tf / (h + Tf);
		A3 = (Kp * ((Tf - Td) * Ti - Tf * Tf) * h) / (Tf * Ti * (h + Tf));
		B3 = (Kp * h) / Ti;
		C3 = (Kp * Td) / Tf;
	}
		
	// Steady state initialization. Assumes e=0.
	// u is the steady state output value.
	void SteadyStateInit(double u) {
		D = 0;
		I = u;
	}

	// Reinitaliztion after paramater change. Keeps the filter state.
	// Can also be used for tracking the states when control is disabled.
	// e is the last value supplied to Update and u the last value returned from Update.
	void ReInit(double e, double u) {
		I = u - C3 * e - D;
	}

	// Full initialization with two sample points.
	void Init(double e0, double u0, double e1, double u1) {
		D = (A1 * u1 + (-A1 * C3 - A1 * B3 - A3) * e1 - A1 * u0 + A1 * C3 * e0) / (A1 - 1);
		I = -(u1 + (-C3 - A1 * B3 - A3) * e1 - A1 * u0 + A1 * C3 * e0) / (A1 - 1);
	}
};

// PID transfer function realization with bilinear integration.
class pid_bil {
protected:
	double e0;
	double D;        // lowpass filter
	double I;        // integrator
	double A1, A3, B3, C3;
public:
	// Update the pid with a new input value 'e'.
	// It advances the pid state by the time step 'h'.
	// When used in real time this function must be called every time 'h' has elapsed.
	// The return value is the control output 'u' of the transfer function.
	double Update(double e) {
		D = A3 * (e + e0) + A1 * D;
		I = B3 * (e + e0) + I;
		e0 = e;
		return C3 * e + I + D;
	};

	// Setup parameters for a serial pid with the transfer function
	// u(s) / e(s) = G * ((Ti * s + 1) / (Ti * s) * (Td * s + 1) * (1 / (Tf * s + 1)));
	// h is the sampling period. It has the same unit as Ti, Td and Tf.
	// G is the pid gain.
	void SerialPid(double h, double G, double Ti, double Td, double Tf) {
		A1 = -(h - 2 * Tf) / (h + 2 * Tf);
		A3 = (G * (Tf - Td) * (Ti - Tf) * h) / (Tf * Ti * (h + 2 * Tf));
		B3 = (G * h) / (2 * Ti);
		C3 = (G * Td) / Tf;
	}

	// Setup parameters for a parallel pid with the transfer function
	// u(s) / e(s) = (Kp + Ki / s + Kd * s) * 1 / (Tf * s + 1)
	// h is the sampling period.
	// Tf is the filter time constant.
	void ParallelPid(double h, double Kp, double Ki, double Kd, double Tf) {
		A1 = -(h - 2 * Tf) / (h + 2 * Tf);
		A3 = -(((Ki * Tf - Kp) * Tf + Kd) * h) / (Tf * (h + 2 * Tf));
		B3 = (Ki * h) / 2;
		C3 = Kd / Tf;
	}

	// Setup parameters for a standard pid with the transfer function
	// u(s) / e(s) = Kp * (1 + 1 / (Ti * s) + s * Td / ((Td / N) * s + 1))
	// h is the sampling period. It has the same unit as Ti and Td.
	// Td/N is the filter time constant.
	// Kp is the pid gain.
	void NStandardPid(double h, double Kp, double Ti, double Td, double N) {
		A1 = -(N * h - 2 * Td) / (N * h + 2 * Td);
		A3 = -(Kp * N * N * h) / (N * h + 2 * Td);
		B3 = (Kp * h) / (2 * Ti);
		C3 = Kp * N + Kp;
	}

	// Setup parameters for a standard pid with the transfer function
	// u(s) / e(s) = Kp * (1 + 1 / (Ti * s) + s * Td) * (1 / (Tf * s + 1))
	// h is the sampling period. It has the same unit as Ti, Td and Tf.
	// Tf is the filter time constant.
	// Kp is the pid gain.
	void StandardPid(double h, double Kp, double Ti, double Td, double Tf) {
		A1 = -(h - 2 * Tf) / (h + 2 * Tf);
		A3 = (Kp * ((Tf - Td) * Ti - Tf * Tf) * h) / (Tf * Ti * (h + 2 * Tf));
		B3 = (Kp * h) / (2 * Ti);
		C3 = (Kp * Td) / Tf;
	}

	// Steady state initialization. Assumes e=0.
	// u is the steady state output value.
	void SteadyStateInit(double u) {
		e0 = 0;
		D = 0;
		I = u;
	}

	// Reinitaliztion after paramater change. Keeps the filter state.
	// Can also be used for tracking the states when control is disabled.
	// e is the last value supplied to Update and u the last value returned from Update.
	void ReInit(double e, double u) {
		I = u - C3 * e - D;
	}

	// Full initialization with two sample points.
	void Init(double e0, double u0, double e1, double u1) {
		this->e0 = e0;
		D = (A1 * u1 + (-A1 * C3 - A1 * B3 - A3) * e1 - A1 * u0 + (A1 * C3 - A1 * B3 - A3) * e0) / (A1 - 1);
		I = -(u1 + (-C3 - A1 * B3 - A3) * e1 - A1 * u0 + (A1 * C3 - A1 * B3 - A3) * e0) / (A1 - 1);
	}
};

template<typename T>
class backcalculation_t : public T {
	double u_min, u_max;
	double cW;
public:
	// Anti windup backcalculation parameters.
	// (min, max) is the allowed output value range.
	// Recommended values for Tw are between
	// [Ti, sqrt(Ti * Td)]
	// If you want simple clipping of the integrator set Tw = h;
	void Backcalculation(double min, double max, double Tw, double h) {
		u_min = min;
		u_max = max;
		cW = Tw > h ? h / Tw : 1;
	}

	double Update(double e) {
		double u = T::Update(e);
		if (u > u_max)
		{
			T::I += cW * (u_max - u);
			return u_max;
		}
		else if (u < u_min)
		{
			T::I += cW * (u_min - u);
			return u_min;
		}
		return u;
	}
};

template<typename T>
class clamping_t : public T {
	double u_min, u_max;
public:
	// Anti windup integrator clamping.
	// (min, max) is the allowed output value range.
	void Clamping(double min, double max) {
		u_min = min;
		u_max = max;
	}

	double Update(double e) {
		double I0 = T::I;
		double u = T::Update(e);
		if (u > u_max)
		{
			if (e > 0)
				T::I = I0;
			return u_max;
		}
		else if (u < u_min)
		{
			if (e < 0)
				T::I = I0;
			return u_min;
		}
		return u;
	}
};
