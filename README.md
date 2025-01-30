# figures_2025_yamda_model_PTC

Matlab and Python scripts for generating the figures in the paper "Phase Resetting in the Yamada Model" by Jacob Ngaha, Neil G. R. Broderick, and Bernd Krauskopf.

The code for both of these programs is written in either [Continuation Core and Toolboxes (COCO)](https://sourceforge.net/projects/cocotools/) (in Matlab) or [AUTO-07P](https://www.github.com/auto-07p/auto-07p/) (in Python and Fortran). 

## Content

Each of the `figX_data.m` COCO scripts makes use of functions in two folders:
- `functions`: Contains the Matlab encoding of the vector fields (`functions/fields/`) and of the boundary conditions (`functions/bcs/`).
- `continuation_scripts`: Contains functions which either calculate an initial solution, setup the COCO continuation problem scturcture, or save the relevant data to a Matlab `.mat` file.

The code is structured as follows:

- `plot_mat_files`  
  After a succesful run of the data generating scripts, the relevant figure data will be saved to `plot_mat_files/figX_data.mat`.

- `fig2`
  - `fig2_data.m`: COCO script for generating the initial periodic orbit solution, and the stable manifold of the stationary point $q$. The initial periodic orbit solution data will also be saved to `plot_mat_files/initial_PO.mat`.
  - `plot_fig2a.m`: Plots Fig. 2(a) - the 3D phase portrait in .
  - `plot_fig2b.m`: Plots Fig. 2(b) - the three components of the Yamada model with respect to time.

- `fig3_fig4`
  - `fig3_fig4_data.m`: COCO script for generating the phase reset for a perturbation $A_{\mathrm{p}}=0.1$, and the phase transition curve (PTC).
  - `plot_fig3.m`: Plots Fig. 3 - the four panel figure depicting the phase portraits (a1 and a2) and temporal trace of the intensity (b1 and b2).
  - `plot_fig4.m`: Plots Fig. 4 - the PTC for the given perturbation.

- `fig5`
  - `fig5_data.m`: COCO script for generating the phase reset for a perturbation $A_{\mathrm{p}}=1.5$, and the phase transition curve (PTC).
  - `plot_fig5.m`: Plots Fig. 5 - the four panel figure depicting the phase portraits (a1 and a2) and temporal trace of the intensity (b1 and b2).

- `fig6`

- `fig7`

- `fig8`


## Usage