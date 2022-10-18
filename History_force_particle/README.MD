# History_force_term
Particle in fluid flow can be tracked in lagrangina form .Particle will be considered as point in the flow and will be travcked using particles equation of motion.
For tracking position, the equation ..................will be used
Particles velcity will be calculated using Maxey Riley equation ...............

On the Left hand side we can see mass x  acceleration  and R.H.S is summation fo forces acting on the particle. F_d is for drag(resistance) force acting on the particles 
as it moves through the fluid and F_B is body force acting on the particle. Other notable force term, which accounts for unsteady in boundary layer develoment around the 
particle is accounted using history force term.This is usually ignored by researchers as it has singular integral and requires storing acceleration of particles.

In 2011, Hinsberg.et.al published a paper https://www.sciencedirect.com/science/article/pii/S0021999110006224?via%3Dihub, which provides an efficient way to numerically implement the history force term.

In this code, we will be looking at an implementation of basset force term using python code . Further, investigation of the influnece of the history force term on the velcoity of particles is also presented using an image. This project was done in collaboration with Dr.Pedram Pakseresht, Professor Sourabh V. Apte
