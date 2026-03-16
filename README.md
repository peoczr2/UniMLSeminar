This is a university research project extending the cabalities of the BEAT model introudeced in Ascarza, Eva, and Ayelet Israeli. "Eliminating unintended bias in personalized policies using bias-eliminating adapted trees (BEAT)." Proceedings of the National Academy of Sciences.
The paper and appendix pdf can be found in the project files.
The code source for the original BEAT can be found on GitHub: https://github.com/ayeletis/beat.git
Our extensions and research project proposal can be found in main file.
Currently there is two possible extensions we will try out: C-BEAT and BTGQ.
We want to first try out our models on simulated data and compare them with the original BEAT, and CF-FD, CF-NP and random allocation to 50%.
We will use the same metrics to measure efficience, Imbalance and delta Policy.
We want multiple different scenarios for the simulation data to well simulate the semi-complex real word messy reality. Obviously incrementally increasing the messyness. But its important to test different scenarios: casuality, omitted variables, confounding, correlations, selection bias, non-stationary historical bias, etc. Some situations are already describe in the paper. But new plausable scenarios reflecting real word application has to be tryed out and tested.