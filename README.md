# lambert-bate

A lightweight Lambert's problem solver. The Universal Variables method introduced in "Fundamentals of Astrodynamics" (Bate, Mueller, and White, 1971) is used. This method was found to be the quickest to converge among those investigated in the review article "Review of Lambert's Problem" ([Sangrà and Fantino, 2021](https://arxiv.org/pdf/2104.05283.pdf)).

## Lambert's problem

[Lambert's problem](https://en.wikipedia.org/wiki/Lambert%27s_problem) is to compute the orbit of a satellite given two points in its orbit and the time of flight between them. It has many uses, from planning spacecraft rendezvous to orbit determination. No formula exists for the orbital parameters as a function for these parameters, but a fairly straightforward formulation exists in which the solution can be found using Newton's method.

Here, we solve only for the satellite's velocities at the two known points in the orbit. If other orbital parameters are required, they can be easily computed from this. For example, the orbital energy is readily available from position and velocity at fixed time, which gives a semi-major axis. The orbital angular momentum is also determined, yielding the semi-latus rectum.

## Benchmarking

10,000 calls of `get_velocities` were run for randomly chosen positions distributed in the cube centered on the origin with side length of 1e11 meters for times of flight between zero and 1e8 seconds. The time of execution and time per function call are compared to algorithms implemented in the Python `lamberthub` library (version 0.1) for reference.

Due to the speed advantages of Rust over Python, `lambert-bate` is more than 50 times faster than the Python implementation.

|Algorithm name | Total time (ms) | Time per function call (μs) |
| --- | --- | --- |
| `lambert-bate` (this crate) | 36| 3.6 |
| `lamberthub.izzo2015` (Python) | 2,700 | 270 |
| `lamberthub.gooding1990` (Python) | 3,700 | 370|
| `lamberthub.avanzini2008` (Python) | 14,000 | 1,400 |

## Change log

- **Version 0.1.0** (June 14, 2022): `lambert-bate` released.