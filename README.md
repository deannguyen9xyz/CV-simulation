# CV-simulation
A Python-based simulation of Cyclic Voltammetry. Models the time-dependent concentration profiles and current response at a planar electrode by numerically solving Fick's Second Law of Diffusion under Nernstian equilibrium.

---

## ▶️ How to Run

Run script `generate_CV_2cycle.py`.

---

## 📊 Result and Conclusion

<img width="887" height="633" alt="image" src="https://github.com/user-attachments/assets/c23f5c76-420e-4183-b7f2-307d99eed1ef" />

*The first cycle of CV doesn't close the loop, while the second one does. This is because the concentration is uniform everywhere (no *R* has been created yet), there is no concentration gradient (**tdC/dx} = 0*) in the initial condition. However, *R* had been created during the first cycle scan, and it is still diffusing out into the solution.

---

## 🎯 Purpose of This Project

* Simulation of CV for a reversible electrochemical system, using Nernst Equation and Fick’s Laws of Diffusion interact at an electrode surface to produce the characteristic "duck-shaped" voltammogram..
* Practicing Explicit Finite Difference Method (FDM) to simulate models the concentration gradients of redox species over time and space..
* GitHub portfolio demonstration.


--- 

## 🧑‍💻 Author

Developed by: Vu Bao Chau Nguyen, Ph.D.

Keywords: Cyclic Voltammetry (CV), Explicit Finite Difference Method (FDM).

---
