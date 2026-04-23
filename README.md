# RSP

MATLAB implementation of **RSP (Robust Sensor Placement)** for source signal extraction from noisy measurements.

This repository is based on the paper:

**Fateme Ghayem, Bertrand Rivet, Rodrigo Cabral Farias, Christian Jutten**  
**“Robust Sensor Placement for Signal Extraction”**  
IEEE Transactions on Signal Processing, 2021

Paper: [RSP_TSP_2021.pdf](https://ghayem.github.io/files/RSP_TSP_2021.pdf)

---

## Overview

This repository implements the method proposed in the paper **“Robust Sensor Placement for Signal Extraction”**, which addresses the problem of placing a limited number of sensors in order to recover a source signal from noisy measurements. The method is designed for settings where both the spatially varying sensor gain and the spatially correlated noise are uncertain, and models both of them as realizations of Gaussian processes. The paper proposes a placement criterion based on maximizing the **probability that the output SNR exceeds a given threshold**, rather than optimizing only an average criterion.

More specifically, the paper considers source extraction from measurements recorded at a subset of candidate sensor positions. Because the number of sensors is limited in many applications, the goal is to choose positions that provide the most robust recovery of the source under uncertainty. The proposed method differs from standard sensor placement methods aimed mainly at field interpolation, and also differs from sensor selection settings where measurements from all positions are already available.

The public repository contains two main implementations:
- `RSP_1D_TSP_2021`
- `RSP_2D_TSP_2021`  
---

## Repository contents

- `RSP_1D_TSP_2021/` — MATLAB code for 1D sensor placement experiments
- `RSP_2D_TSP_2021/` — MATLAB code for 2D sensor placement experiments

The repository is primarily MATLAB code, and the two main folders suggest separate implementations for one-dimensional and two-dimensional experimental settings.

---

## Problem setting

The paper studies a scenario where a source signal propagates through a structure and is observed by several sensors. Since the number of sensors is limited, one must decide where to place them among candidate locations. The objective here is specifically **signal extraction**, meaning that the sensor positions should help recover the source as reliably as possible from noisy measurements.

A key aspect of the method is that uncertainty is explicitly modeled:

- the spatial sensor gain is uncertain,
- the noise is spatially correlated and uncertain,
- the SNR is therefore also uncertain.

To account for this, the paper models the gain and noise processes with Gaussian processes and proposes a robust criterion that maximizes the probability that the resulting SNR is above a prescribed threshold.

---

## Main idea of the method

Many sensor placement methods optimize an average quantity, such as a mean SNR or an interpolation-oriented score. In contrast, this paper proposes a **robust probabilistic criterion**:

- instead of maximizing only the expected SNR,
- it maximizes the probability that the SNR is larger than a fixed threshold.

Because the gain and noise are modeled as Gaussian processes, this probability can be computed analytically or efficiently approximated within the proposed framework. The paper also introduces a **sequential maximization strategy**, where sensor positions are chosen one by one, in order to reduce the computational cost of full joint optimization over all sensor locations.

In summary, the method has two main strengths:

1. **Robustness to uncertainty** in gain and correlated noise  
2. **Reduced computational complexity** through sequential sensor selection rather than exhaustive joint search

---

## 1D and 2D implementations

The repository contains two experiment folders:

- `RSP_1D_TSP_2021`
- `RSP_2D_TSP_2021`

These correspond to one-dimensional and two-dimensional spatial sensor placement experiments, matching the types of numerical studies discussed in the paper. The paper presents numerical results demonstrating the method in simulated settings and compares the proposed robust strategy against alternative criteria.

This split is useful for studying the method in different spatial geometries.

---

## Method summary

The proposed approach can be summarized as follows:

1. Model the uncertain spatial sensor gain as a Gaussian process  
2. Model the spatially correlated noise as a Gaussian process  
3. Define the output SNR at candidate sensor positions  
4. Compute the probability that the SNR exceeds a desired threshold  
5. Choose the sensor positions that maximize this probability  
6. Use a sequential strategy to add sensors one at a time when full joint optimization is too expensive

This criterion is designed for **robust source extraction**, not just interpolation accuracy.

---

## Why this method is different

The work distinguishes this problem from two related settings:

- **Sensor selection**: all sensor signals are already available, and one selects a subset afterward
- **Sensor placement**: one must decide where to place the sensors before measurements are available

The method in this repository belongs to the second category. It is also different from criteria that only maximize average SNR, because it accounts for uncertainty and optimizes a probabilistic robustness criterion.

---

## Expected workflow

A typical workflow with this repository is:

1. Open either the `RSP_1D_TSP_2021` or `RSP_2D_TSP_2021` folder
2. Add the folder and its subfolders to the MATLAB path
3. Run the main experiment scripts in that folder
4. Compare the selected sensor positions and the resulting robustness metrics

---

## Experimental results in the paper

The paper reports that the proposed robust sensor placement strategy shows **superior robustness** compared with:
- standard sensor placement criteria aimed at interpolating the spatial gain,
- and a previously proposed criterion that maximizes the average SNR.

The numerical experiments are designed to show that optimizing the probability of exceeding an SNR threshold leads to more reliable source extraction under uncertainty.

---

## Citation

If you use this code, please cite:

```bibtex
@article{ghayem2021rsp,
  title={Robust Sensor Placement for Signal Extraction},
  author={Ghayem, Fateme and Rivet, Bertrand and Cabral Farias, Rodrigo and Jutten, Christian},
  journal={IEEE Transactions on Signal Processing},
  year={2021}
}
