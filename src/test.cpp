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


#define _CRT_SECURE_NO_WARNINGS
#include "cmon-pid.h"
#include <stdio.h>
#include <array>
#include <functional>
#include <fstream>

using namespace std;

// A vector of doubles that supports some arithmetic
template<size_t n>
class darray : public array<double, n> {
	typedef darray<n> _ty;
public:
	darray(const array<double, n>& x) : array<double, n>(x) {}
	darray() {}
	_ty operator*(_ty b) {
		_ty r;
		for (size_t i = 0; i < n; ++i) r[i] = _ty::operator[](i) * b[i];
		return r;
	}
	_ty operator*(double b) {
		_ty r;
		for (size_t i = 0; i < n; ++i) r[i] = _ty::operator[](i) * b;
		return r;
	}
	_ty operator+(_ty b) {
		_ty r;
		for (size_t i = 0; i < n; ++i) r[i] = _ty::operator[](i) + b[i];
		return r;
	}
	_ty operator+(double b) {
		_ty r;
		for (size_t i = 0; i < n; ++i) r[i] = _ty::operator[](i) + b;
		return r;
	}
};

// explicit 4th order runge kutta
template<typename T>
void step_rk4(T& sys, double dt)
{
	auto X = sys.X;
	auto dX1 = sys.F() * dt;
	sys.t = sys.t + dt  * 0.5;
	sys.X = X + dX1 * 0.5;
	auto dX2 = sys.F() * dt;
	sys.X = X + dX2 * 0.5;
	auto dX3 = sys.F() * dt;
	sys.t = sys.t + dt * 0.5;
	sys.X = X + dX3;
	auto dX4 = sys.F() * dt;
	sys.X = X + (dX1 + dX2 * 2.0 + dX3 * 2.0 + dX4) * (1 / 6.0);
}

// A linear continuous transfer function system of order n in controlable canonical form
template<size_t n>
class tf {
public:
	typedef darray<n> state_t;
	typedef darray<n> state_diff_t;
	typedef darray<n + 1> coef_t;
	coef_t a;
	coef_t b;

	state_t X;  // current state
	double t;   // current time, but the system is time invariant.
	double u;	// input variable

	// constructor
	// a0 is the denominator polynom of the transfer function
	// b0 is the numerator polynom of the transfer function
	tf(array<double, n + 1> a0, array<double, n + 1> b0) : a(a0), b(b0), u(0), t(0) 
	{
		static_assert(n <= 6,"The polynomial coefficients used by the controlable canonical form are ill-conditioned "
			"and may not work for high-order systems.");
		b = b * (1 / a[n]);
		a = a * (1 / a[n]); // Make the denominator polynom monic.
		X.fill(0);
	}

	void SteadyStateInit(double y) {
		t = 0;
		X.fill(0);
		X[0] = y / b[0];
		u = a[0] * y / b[0];
	}

	double Xn() {
		double xn = u;
		for (auto i = 0; i < n; ++i)
			xn += -a[i] * X[i];
		return xn;
	}

	// returns the current 1st derivative of the state variables.
	state_diff_t F()
	{
		state_diff_t dy;
		for (auto i = 0; i < n - 1; ++i)
			dy[i] = X[i + 1];
		dy[n - 1] = Xn();
		return dy;
	}

	// returns the current output value of the transfer function.
	double Y() {
		double y = b[n] * Xn();
		for (auto i = 0; i < n; ++i) {
			y += b[i] * X[i];
		}
		return y;
	}
};


class simfile {
	bool isEol;
	ofstream o;
public:
	void eol() {
		o << '\n';
		isEol = true;
	}
	void sep() {
		if (isEol)
			isEol = false;
		else
			o << ',';
	}
	void open(const char* filename)
	{
		o.precision(18);
		isEol = true;
		if (o.is_open())
			o.close();
		o.open(filename, o.out | o.trunc);
	}
	simfile& operator<<(const string& s) { sep();  o << s; return *this; }
	simfile& operator<<(double d) { sep();  o << d;  return *this; }
};

struct simulationItem {
	struct Con {
		double setpoint = 0;
		double distortion = 0;
	};
	virtual void Title(simfile& out) = 0;
	virtual void Init(Con& con, simfile& out) = 0;
	virtual void Step(Con& con, simfile& out, double dt) = 0;
};

template<typename... Arguments>
void simulate(double duration, double dt, const char* filename, Arguments&... args)
{
	constexpr auto n{ sizeof...(Arguments) };
	simulationItem* isim[] = { &args... };

	simfile out;
	double t = 0;
	simulationItem::Con con;

	out.open(filename);
	out << "t";
	for (auto sim : isim)
		sim->Title(out);
	out.eol();

	out << t;
	for (auto sim : isim)
	{
		sim->Init(con, out);
	}
	out.eol();
	while (t < duration)
	{
		t = t + dt;
		out << t;
		for (auto sim : isim)
		{
			sim->Step(con , out, dt);
		}
		out.eol();
	}
}

// Simulates a step function for the setpoint
class StepSim : public simulationItem {
	double t;
	bool bDisturb;
public:
	StepSim(bool b) : t(0), bDisturb(b) {}
	void Title(simfile& out) override {
		out << "setpoint";
	}

	void Init(Con& con, simfile& out) override {
		t = 0;
		con.setpoint = 0;
		con.distortion = 0;
		out << con.setpoint;
	}

	void Step(Con& con, simfile& out, double dt) override
	{
		t = t + dt;
		if (t >= 0.5)
			if (bDisturb)
				con.distortion = 1;
			else
				con.setpoint = 1;
		out << con.setpoint;
	}
};

// Simulates a transfer function system
template<typename TF>
class TfSim : public simulationItem
{
	TF tf;
	const char* hdr;

public:
	TfSim(const TF& tf0, const char* hdr0) : tf(tf0), hdr(hdr0) {}
	void Title(simfile& out) override {
		out << hdr;
	}

	void Init(Con& con, simfile& out) override {
		tf.SteadyStateInit(con.setpoint);
		out << con.setpoint;
	}

	void Step(Con& con, simfile& out, double dt) override
	{
		tf.u = con.setpoint;
		step_rk4(tf, dt);
		out << tf.Y();
	}
};

// Simulates a cmon-pid
template<typename PID, typename TF>
class PidLoopSim : public simulationItem
{
	PID pid;
	TF tf;
	double h;
	double clock;
	string hdr;

public:
	PidLoopSim(const PID& p0, const TF& tf0, double h0, string hdr0) :
		hdr(hdr0), pid(p0), tf(tf0), h(h0), clock(0) {}

	void Title(simfile& out) override {
		out << "u" + hdr;
		out << "y" + hdr;
	}

	void Init(Con& con, simfile& out) override {
		tf.SteadyStateInit(con.setpoint);
		pid.SteadyStateInit(tf.u);
		out << tf.u;
		out << con.setpoint;
	}

	void Step(Con& con, simfile& out, double dt) override
	{
		step_rk4(tf, dt);
		double y = tf.Y();

		clock += dt;
		if (clock >= h)
		{
			tf.u = pid.Update(con.setpoint - y) + con.distortion;
			clock -= h;
		}
		out << tf.u;
		out << y;
	}
};

template<typename T> void antiwindup(backcalculation_t<T>& p, double min, double max, double Tw, double h)
{
	p.Backcalculation(min, max, Tw, h);
};

template<typename T> void antiwindup(clamping_t<T>& p, double min, double max, double Tw, double h)
{
	p.Clamping(min, max);
};

struct SimulationParameters {
	double dt = 1.0 / 1000;	// simulation step
	double duration = 10;
	double ymax = 20000000;
	bool bDisturbance = false;

	double T1 = 8;        // First time constant of test system
	double T2 = 3;        // Second time constanf of test system
	double k = 1;         // DC gain of test system
	double Tc = 1;        // Rise time objective of closed loop

	// pid parameters for the serial pid (G*(Td*s+1)*(Ti*s+1))/(Ti*s*(Tf*s+1))
	double G = T1 / k / Tc;  
	double Ti = T1;
	double Td = T2;

	double Tf = 0.1;
	double h = 0.1;
	double Tw = h;

	template< template<typename T> class Antiwindup >
	void sim(const char* filename)
	{
		// Test system
		// Transfer function of a second order system with time constants at T1 and T2 and gain k. k/((T1*s+1)*(T2*s+1))
		tf<2> openloop({ 1, T2 + T1, T1 * T2 }, { k, 0, 0 });

		// closed loop transfer function of the serial pid (G*(Td*s+1)*(Ti*s+1))/(Ti*s*(Tf*s+1))
		// with the second order test system k/((T1*s+1)*(T2*s+1))
		tf<4> closedloop(
			{ G * k, (G * Ti + G * Td) * k + Ti, G * Td * Ti * k + (Tf + T2 + T1) * Ti, ((T2 + T1) * Tf + T1 * T2) * Ti, T1 * T2 * Tf * Ti },
			{ G * k, (G * Ti + G * Td) * k, G * Td * Ti * k }
		);

		// PID backward euler
		Antiwindup<pid_bwe> pid1;
		pid1.SerialPid(h, G, Ti, Td, Tf);
		antiwindup(pid1, -ymax, ymax, Tw, h);

		// PID tustin
		Antiwindup<pid_bil> pid2;
		pid2.SerialPid(h, G, Ti, Td, Tf);
		antiwindup(pid2, -ymax, ymax, T2, h);

		auto s0 = StepSim(bDisturbance);
		auto s1 = TfSim(closedloop, "continuous");
		auto s2 = TfSim(openloop, "openloop");
		auto p1 = PidLoopSim(pid1, openloop, h, "_bwe");
		auto p2 = PidLoopSim(pid2, openloop, h, "_bil");
		simulate(duration, dt, filename, s0, s1, s2, p1, p2);
	}

};

int main()
{
	int umax = 10;
	{
		SimulationParameters s;
		s.Tf = 0.1;
		s.sim<clamping_t>("sim_fast_sample.csv");
	}

	{
		SimulationParameters s;
		s.Tf = 0.1;
		s.h = 0.5;
		s.sim<clamping_t>("sim_slow_sample.csv");
	}

	{
		SimulationParameters s;
		s.Tf = 0.5;
		s.sim<clamping_t>("sim_low_bandwidth.csv");
	}

	{
		SimulationParameters s;
		s.Tf = s.h / 100;
		s.sim<clamping_t>("sim_high_bandwidth.csv");
	}

	{
		SimulationParameters s;
		s.duration = 50;
		s.bDisturbance = true;
		s.sim<clamping_t>("sim_disturbance_rejection.csv");
	}

	{
		// Backcalculation with Tw = Ti
		SimulationParameters s;
		s.duration = 50;
		s.ymax = umax;
		s.Tw = s.Ti;  // back calculation factor
		s.sim<backcalculation_t>("sim_backcalc_ti.csv");
	}

	{
		// Integrator clamping
		SimulationParameters s;
		s.duration = 50;
		s.ymax = umax;
		s.sim<clamping_t>("sim_clamping.csv");
	}

	{
		// Improved disturbance rejection. Step response.
		SimulationParameters s;
		s.Ti = 2;
		s.sim<clamping_t>("sim_Ti2.csv");
	}

	{
		// Improved disturbance rejection.
		SimulationParameters s;
		s.Ti = 2;
		s.duration = 50;
		s.bDisturbance = true;
		s.sim<clamping_t>("sim_Ti2_disturbance.csv");
	}

	{
		// Clamping with improved disturbance recejction.
		SimulationParameters s;
		s.duration = 50;
		s.Ti = 2;
		s.ymax = umax;
		s.sim<clamping_t>("sim_Ti2_clamping.csv");
	}

	return 0;
}

