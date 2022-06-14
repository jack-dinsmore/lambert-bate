# lambert-bate

A Lambert's problem solver using the Universal Variables method introduced in "Fundamentals of Astrodynamics" (Bate, Mueller, and White, 1971). This algorithm was found to be the quickest to converge among those investigated in the review article "Review of Lambert's Problem" ([Sangrà and Fantino, 2021](https://arxiv.org/pdf/2104.05283.pdf))

## Lambert's problem

[Lambert's problem](https://en.wikipedia.org/wiki/Lambert%27s_problem) is to compute the orbit of a satellite given two points in its orbit and the time of flight between them. It has many uses, from planning spacecraft rendezvous to orbit determination. No formula exists for the orbital parameters as a function for these parameters, but a fairly straightforward numerical algorithm to extract them does exist.

Here, we solve only for the satellite's velocities at the two known points in the orbit. If other orbital parameters are required, they can be easily computed from this. For example, the orbital energy is readily available from position and velocity at fixed time, which gives a semi-major axis. The orbital angular momentum is also determined, yielding the semi-latus rectum.