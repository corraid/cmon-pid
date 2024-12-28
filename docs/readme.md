
# cmon-pid
A implementation of several common PID controller transfer functions in C++, header only, without dependencies.

* Realized PID transfer functions
  + [Serial](index.html#SeriesPID)
  + [Parallel](index.html#ParallelPID)
  + [N-Standard](index.html#NStdPID)
  + [Standard](index.html#StdPID)

* [Numerical solution](index.html#TFSolutions)
  + Backward euler
  + Bilinear / Tustin
  
* Anti-windup
  + [Clamping](index.html#SimClamping) 
  + [Back-calculation](index.html#SimBackcalc)
  
* Initialization
  + Steady State
  + Bump-less parameter change
  + Sampling points

* [Simulation](index.html#Simulation) of a pid in closed loop with a second order test system

Usage:
~~~
#include <cmon-pid.h>
extern double GetActuator();
extern double Setpoint();
extern double SensorReading();
extern void SetActuator(double);
extern void Sleep(double);

void PidTask()
{

	// Select a pid object of class pid_bwe (backward euler)
	// or pid_bil (tustin) and
	// a anti-windup template class backcalculation_t
	// or clamping_t
	clamping_t<pid_bwe> pid;

	// Set PID parameters
	constexpr double samplingTime = 0.1;
	constexpr double gain = 10;
	constexpr double T1 = 10;
	constexpr double T2 = 3;
	pid.SerialPid(samplingTime, gain, T1, T2, samplingTime / 2);

	// Set clamping range
	pid.Clamping(-2, 10);

	// Initialize it
	pid.SteadyStateInit(GetActuator());

	// Control loop
	for (;;)
	{
		double e = Setpoint() - SensorReading();
		double y = pid.Update(e);
		SetActuator(y);
		Sleep(samplingTime);
	}
}
~~~
