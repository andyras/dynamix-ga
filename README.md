## Creating a new set of variables to optimize
In order to optimize a new set of variables, these things need to be added to the source code:

1. `initializer.cpp` needs a new function to take an individual's genome and change the proper parameters in the `dynamix` *Params* object. The naming convention for the function is `init_XXX`, where `XXX` matches the `initializer` parameter from the `dynamix-ga` *GAParams* object.

2. `main.cpp` needs to have an `else if` condition added where the proper upper and lower bounds for each parameter are set according to the value of `GAParams.initializer`.

3. `objective.cpp` needs to be updated in the `singleObjective` and `doubleObjective` functions to call the proper `init_XXX` function according to the value of `GAParams.initializer`.

4. `README.md` should be updated!

I would like to make a more flexible way to change parameters, but a lot of that depends on how things are initialized in the original `dynamix` code... some things (like changing the number of states on the acceptor) are a bit tricky to set up.


## Input parameters

The file `ga.in` contains input parameters relevant for the genetic algorithm.


### `objectiveType`
This parameter whether the objective function is `single` (obtains an observable from one run of `dynamix`) or `double` (obtains an observable from two runs of `dynamix`).
Values:

* `single`

* `double`


### `objective`
This parameter determines what observable is computed by the objective function.
Values:

* `acceptorPeak`: returns the peak value of electron population on acceptor

* `acceptorAvg`: returns average value of electron population on acceptor

* `acceptorAvgAfterPeak`: returns average of electron population on acceptor after peak value of population is reached

* `acceptorFinal`: returns electron population on acceptor at the end of the simulation


### `doubleObjectiveType`
This parameter determines the difference between the two runs of `dynamix` that are compared when using a double objective function.
Values:

* `coherence`: the first propagation is for a coherent starting state, and the second is incoherent (only diagonal elements of density matrix are populated to start)


### `initializer`
This parameter determines the variables which are subject to obtimization by the GA.

Values:

* `torsion`
The static component (`torsionSin2V0`) and oscillating component (`torsionSin2V1`) of the torsionally-modulated coupling, and the bridge sites' energy (both are assumed to be the same energy) are varied.

* `g1g2g1_c`
The relaxation rate on the donor (`gamma1`) and acceptor (`gamma1_c`), and the decoherence rate (`gamma2`) are varied.

* `wavepacket`
The central energy (`bulkGaussMu`) and width (`bulkGaussSigma`) of the initial electronic wavepacket on the donor are varied.

* `wavepacketGammas`
Combined variables from `g1g2g1_c` and `wavepacket` options.


### `minmax`
This parameter determines whether the GA minimizes or maximizes the objective function.

Values:

* `min`

* `max`


### `popsize`
This parameter determines the size of the population subject to evolution.

Values: int > 1


### `pMut`
This parameter determines the probability of mutation of each gene when creating a new individual.

Values: float [0,1]


### `pCross`
This parameter determines the crossover probability when creating a new individual.

Values: float [0,1]


### `convergence`
This parameter determines the convergence criterion for the relative change in the objective function of the best individual in each generation. When the change is below this threshold for 20 generations the GA terminates.