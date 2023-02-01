# seasize
A comprehensive numerical sizing script in MATLAB for motor, gearbox (and spring) for (series elactic) actuators.

Designing an actuator is challenging. Often, only nominal torque and speed requirements are met, and acceleration requirements and inertia are ignored, which can result in an undersized actuator, unable to perform the required task.
This project intends to help actuator designers and researchers properly size an actuator for a well-defined task.
The code implements a sizing method introduced by van de Straete et al. [2], and expands upon it by adding the SEA's spring to the evaluation.

## Example: Combined reduction ratio and stiffness evaluation
The matlab script file "[example.m](https://github.com/sadegroo/seasize/blob/main/oo_scripts/example.m)" contains an example use of the [classes](https://github.com/sadegroo/seasize/tree/main/classes) capable of evaluating the effect of reduction and stiffness simultaneously.

### Example output
#### Motion profile
![Motion profile](https://github.com/sadegroo/seasize/blob/main/figs/gaitprofile31.png)
#### Reduction ratio evaluation
![Reduction ratio evaluation](https://github.com/sadegroo/seasize/blob/main/figs/geareval31.png)
#### Relative RMS torque with spring
![Relative RMS torque with spring](https://github.com/sadegroo/seasize/blob/main/figs/relrmstor31.png)
#### Relative peak absolute power with spring
![Relative peak absolute power with spring](https://github.com/sadegroo/seasize/blob/main/figs/relpkabspow31.png)
#### Relative average power with spring (proportional to energy use)
![Relative average power with spring](https://github.com/sadegroo/seasize/blob/main/figs/relavgpow31.png)

## Code used in paper
This project contains legacy code, as it was used at the time we designed our SEA and wrote the paper. It is containded in the [mainscripts](https://github.com/sadegroo/seasize/tree/main/mainscripts) folder. This code does not evaluate a N-K grid (reduction-stiffness), but rather evalutes the reduction ratio first. Afterwards, a ratio needs to be selected, and the effect of stiffness is evaluated. Finally, the motor limits need to be checked. This procedure can be iterative and was improved upon by using classes and the grid approach.

## Resources
- [1] [<div class="csl-entry">De Groof, S., Zhang, Y., Peyrodie, L., Sevit, R., vander Poorten, E., Aertbelien, E., &#38; Labey, L. (2022). Design of a Series Elastic Actuator for an Assistive Exoskeleton using Numerical Methods and Gait Data. <i>ACTUATOR 2022; International Conference and Exhibition on New Actuator Systems and Applications</i>, 1–4.</div>](https://ieeexplore-ieee-org.kuleuven.e-bronnen.be/document/9899227)
- [2] [<div class="csl-entry">van de Straete, H. J., Degezelle, P., Schutter, J. de, &#38; Belmans, R. J. M. (1998). Servo motor selection criterion for mechatronic applications. <i>IEEE/ASME Transactions on Mechatronics</i>, <i>3</i>(1), 43–50. https://doi.org/10.1109/3516.662867</div>](https://ieeexplore-ieee-org.kuleuven.e-bronnen.be/document/662867)


