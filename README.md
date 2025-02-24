# figures_2025_yamda_model_PTC

Matlab and Python scripts for generating the figures in the paper "Phase Resetting in the Yamada Model" by Jacob Ngaha, Neil G. R. Broderick, and Bernd Krauskopf.

The code for both of these programs is written in either [Continuation Core and Toolboxes (COCO)](https://sourceforge.net/projects/cocotools/) (in Matlab) or [AUTO-07P](https://www.github.com/auto-07p/auto-07p/) (in Python and Fortran).

## Usage

To generate the data and plot the figures from this paper, make sure you have the latest version of both COCO and AUTO installed.

Simply run the COCO and `plot_figX.m` scripts in Matlab. The AUTO scripts can be run with:
```sh
auto figX_data.py
```
You can also run AUTO through the `auto` command,
```sh
auto figX_data.py
```
or, if you have configured AUTO as a python package, with the standard python:
```sh
python figX_data.py
```
To do this through Conda, just run:
```sh
conda develop /path/to/auto/07p/python
conda develop /path/to/auto/07p/python/auto
```
There are two lines at the top of each AUTO code which can also add the auto commands to the Python Path:
```python
import sys
sys.path.append('/path/to/auto/07p/python')
sys.path.append('/path/to/auto/07p/python/auto')
```
Just change `/path/to/` to the directory where you have AUTO installed.

## Files

Each of the `calc_figX_data.m` COCO scripts makes use of functions in two folders:
- `functions`: Contains the Matlab encoding of the vector fields (`functions/fields/`) and of the boundary conditions (`functions/bcs/`).
- `continuation_scripts`: Contains functions which either calculate an initial solution, setup the COCO continuation problem scturcture, or save the relevant data to a Matlab `.mat` file.

In the main directory, there is also a Matlab function file `set_figure_dimensions.m`. This function defines the dimensions of the axis 'Position' property to a specific dimension in centimetres. For the paper, I import the figures without any labels into Inkscape, and use the textbox tool in Inkscape to generate all of the labels, ticks, and text. `set_figure_dimensions` therefore allows me much finer control over the size of the figures.

This repository is structured as follows:

- `data_files`  
  After a succesful run of the data generating scripts, the relevant figure data will be saved to `data_files/figX_data.mat`.

- `fig2`
  - `calc_fig2_data.m`: COCO script for generating the initial periodic orbit solution, and the stable manifold of the stationary point $q$.
  - `plot_fig2a.m`: Plots Fig. 2(a) - the 3D phase portrait in .
  - `plot_fig2b.m`: Plots Fig. 2(b) - the three components of the Yamada model with respect to time.
  - `save_fig2_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig3`
  - `calc_fig3_data.m`: COCO script for generating the phase reset for a gain perturbation -  $\mathbf{\mathit{d}} = (1, 0, 0)$ - with $A_{\mathrm{p}} = 0.1$, and the phase transition curve (PTC).
  - `plot_fig3a1.m`: Plots Fig. 3(a1) - the phase portrait of the perturbed orbit in the $G-I$ for $A_{\mathrm{p}} = 0.1$ at $\vartheta_{\mathrm{o}} = 0.0$.
  - `plot_fig3a2.m`: Plots Fig. 3(a2) - the phase portrait of the perturbed orbit in the $G-I$ for $A_{\mathrm{p}} = 0.1$ at $\vartheta_{\mathrm{o}} = 0.3$.
  - `plot_fig3b1.m`: Plots Fig. 3(b1) - Temporal trace of the intensity for the phase reset in Fig. 3(a1).
  - `plot_fig3b2.m`: Plots Fig. 3(b2) - Temporal trace of the intensity for the phase reset in Fig. 3(a2).
  - `save_fig3_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig4`
  - `calc_fig4_data.m`: COCO script for generating the phase reset for a gain perturbation -  $\mathbf{\mathit{d}} = (1, 0, 0)$ - with $A_{\mathrm{p}} = 1.5$, and the phase transition curve (PTC).
  - `plot_fig4a1.m`: Plots Fig. 4(a1) - the phase portrait of the perturbed orbit in the $G-I$ for $A_{\mathrm{p}} = 1.5$ at $\vartheta_{\mathrm{o}} = 0.0$.
  - `plot_fig4a2.m`: Plots Fig. 4(a2) - the phase portrait of the perturbed orbit in the $G-I$ for $A_{\mathrm{p}} = 1.5$ at $\vartheta_{\mathrm{o}} = 0.3$.
  - `plot_fig4b1.m`: Plots Fig. 4(b1) - Temporal trace of the intensity for the phase reset in Fig. 5(a1).
  - `plot_fig4b2.m`: Plots Fig. 4(b2) - Temporal trace of the intensity for the phase reset in Fig. 5(a2).
  - `save_fig4_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig5`
  - `calc_fig5_data.m`: COCO script for generating the PTCs for gain perturbations - $\mathbf{\mathit{d}} = (1, 0, 0)$ - for two amplitudes: $A_{\mathrm{p}} = 0.1$ and $1.5$.
  - `plot_fig5.m`: Plots Fig. 5 - PTC for two perturbation amplitudes.
  - `save_fig5_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig6`
  - `calc_fig6_data.m`: COCO script for generating the PTCs for a few different perturbation amplitudes: $A_{\mathrm{p}} = 0.05, 0.1, 0.2, 0.55, 1.0, 1.5$, and $2.0$.
  - `plot_fig6a1.m`: Plots Fig. 6(a1) - the phase portrait of the shifted periodic orbits for the "weaker" perturbation ampltiude, $A_{\mathrm{p}} \leq A_{\mathrm{p}}^{*} \approx 0.55$
  - `plot_fig6a2.m`: Plots Fig. 6(a2) - the phase portrait of the shifted periodic orbits for the "stronger" perturbation ampltiude, $A_{\mathrm{p}} \geq A_{\mathrm{p}}^{*} \approx 0.55$
  - `plot_fig6b1.m`: Plots Fig. 6(b1) - PTCs of the "weaker" perturbation ampltiude, $A_{\mathrm{p}} \leq A_{\mathrm{p}}^{*} \approx 0.55$
  - `plot_fig6b2.m`: Plots Fig. 6(b2) - PTCs of the "stronger" perturbation ampltiude, $A_{\mathrm{p}} \geq A_{\mathrm{p}}^{*} \approx 0.55$
  - `save_fig6_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig7`
  - `calc_fig7_data.py`: AUTO script for generating a whole lot of PTCs for increasing perturbation amplitudes for a $G$-perturbation.
  - `plot_fig7a.m`: Plots Fig. 7(a) - Surface of PTCs from one viewpoint.
  - `plot_fig7b.m`: Plots Fig. 7(b) - The same surface of PTCs from another viewpoint.
  - `save_fig7_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig8`
  - `calc_fig8_data.m`: COCO script for generating the phase reset for an intensity perturbation - $\mathbf{\mathit{d}} = (0, 0, 1)$ - with $A_{\mathrm{p}} = 1.5$, and the phase transition curve (PTC).
  - `plot_fig8a1.m`: Plots Fig. 8(a1) - the phase portrait of the perturbed orbit in the $G-I$ for $A_{\mathrm{p}} = 1.5$ at $\vartheta_{\mathrm{o}} = 0.0$.
  - `plot_fig8a2.m`: Plots Fig. 8(a2) - the phase portrait of the perturbed orbit in the $G-I$ for $A_{\mathrm{p}} = 1.5$ at $\vartheta_{\mathrm{o}} = 0.3$.
  - `plot_fig8b1.m`: Plots Fig. 8(b1) - Temporal trace of the intensity for the phase reset in Fig. 8(a1).
  - `plot_fig8b2.m`: Plots Fig. 8(b2) - Temporal trace of the intensity for the phase reset in Fig. 8(a2).
  - `save_fig5_data.m`: Saves the figure data to a Matlab `.mat` file.

- `fig9`
  - `calc_fig9_data.py`: AUTO script for generating a whole lot of PTCs for increasing perturbation amplitudes for an $I$-perturbation.
  - `plot_fig9.m`: Plots Fig. 9 - Surface of PTCs from one viewpoint.
  - `save_fig9_data.m`: Saves the figure data to a Matlab `.mat` file.

- `set_figure_dimensions.m`: Matlab function file to define the dimensions of the 'Position' property of the figure axis in a way that works for me :)