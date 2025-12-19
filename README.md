# Axial_Groove_Heat_Pipe_Capillary_Limit

**Adiabatic-section Flow in Axial-Groove Heat Pipes with Varying Meniscus Curvature**

Code Repository for the Associated Publication

**Authors:**  
Marc Hodes¹*, Haotian Jia¹, Andrew Daetz¹, and Toby Kirk²  
¹ Department of Mechanical Engineering, Tufts University, Medford, MA 02155, USA  
² School of Mathematical Sciences, University of Southampton, Southampton, SO17 1BJ, UK  

\*Corresponding author: marc.hodes@tufts.edu

---

## Overview

This repository contains the MATLAB code used to generate the results presented in the paper for the limiting case when the protrusion angle defined in the referenced journal publication is zero degrees.

The main code file is:

- **Order_0th_Solving_v68_Referee.m**

To use the code, please refer to the paper and review the comments provided within the scripts.

---
## Code Structure
The script is modularly organized into clearly marked sections, enabling reproduction of individual plots.

Main Script: The primary file follows the logical flow of the analysis. Each major computational step—solving for $\tilde{\xi}$ and $\tilde{Q}_v$, determining the capillary limit, calculating $p^$ vs $z^$, etc. Sections are explicitly named at their beginning.

Functions (FUNCTION SECTION): Reusable code is organized into functions. Some are defined at the end of the main script within a dedicated FUNCTION SECTION, while others are housed in separate .m files. Files beginning with **“A_”** contain analytical functions that are called by the main script. All serve the same purpose: to modularize operations and keep the main script clean.

## How to Reproduce a Plot
To generate a specific figure:

Set Fluid Properties & Solver Parameters: In the first section, configure settings such as temperature (e.g., 60, 200, or 300). Then define the Geometry inn the second section, specify geometric parameters (e.g., set $R_g^*$ to 1.18, 1.20, or 1.22 mm).

Execute Sections Sequentially: Follow the in-line comments and run the code section-by-section in the provided order. Most sections depend on outputs from previous ones (e.g., the capillary limit calculation requires the solutions for $\tilde{\xi}$ and $\tilde{Q}_v$).

By proceeding stepwise through these defined sections, users can seamlessly reconstruct the analysis pipeline and reproduce any plot from the study.

---

## Contact

For questions about the code or methodology, please contact:

- **Marc Hodes** – marc.hodes@tufts.edu  
- **Haotian Jia** – haotian.jia@tufts.edu

---

## Citation

If you use this code in academic work, please cite the associated publication.

