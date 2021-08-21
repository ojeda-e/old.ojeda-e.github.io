---
title:  "GSoC Report: Curvature for MDA Universes"
categories:
    - Blog
tags:
    - Curvature
    - Biological Membranes
    - Lipid Bilayers
toc: true
---


A summary of my GSoC journey 

------------------------------
The main goal of my [GSoC project](https://summerofcode.withgoogle.com/projects/#5098282306502656)
was to develop a tool to calculate membrane curvature
from MD simulations using [MDAnalysis]. We
aimed to have the following features in the membrane curvature tool:

- Surfaces derived from an `AtomGroup` of reference.
- Calculation of mean and Gaussian curvature.
- Implement multiframe and averaged-over-frames analysis.
- Plug-and-play with visualization libraries to obtain 2D curvature profiles.
- Data visualization made easy.

![patches](../../assets/images/patches.png)
 

# Why Membrane Curvature?

In the wide range of tools that are available to analyze Molecular Dynamics (MD)
simulations, user-friendly, actively-maintained, and well-documented tools to
calculate membrane curvature are still difficult to find. I was motivated to
share a tool to calculate membrane curvature that I initially developed as part
of my PhD at the [Biocomputing Group](https://ucalgary.ca/biocomputing/home).
Membrane curvature is a phenomenon that can be investigated via MD simulations, and
I had an interest to share this tool with the wider MD community. 


# Contributions
Keeping in mind the goals and motivations behind my GSoC project,
I would like to hightlight three areas of contributions that were key to the
membrane-curvature MDAnalysis tool: [Core functions](#core-functions),
[AnalysisBase](#Analysis-base), and [Documentation](#documentation).

## Core functions ([#34], [#40], [#44])
The initial version of the code ([#9]) contained functions to map elements from the `AtomGroup` 
of reference into a grid of dimensions defined by the simulation box. [#9] also included the functions to
calculate mean and Gaussian curvature were also included. After refactoring ([#34]), the core functions 
of MembraneCurvature were tuned-up and cleaned.

- Derive surface and normalize grid ([#40])
- Calculate mean and Gaussian curvature. ([#44])

These functions are our bricks to build the membrane-curvature AnalysisBase.


## Analysis Base ([#43], [#48])

In the
[MembraneCurvatureAnalysisBase](https://github.com/MDAnalysis/membrane-curvature/blob/main/membrane_curvature/base.py),
we define the initial arrays for surface, mean, and Gaussian curvature in the
`_prepare()` subclass.  In `_single_frame()`, AnalysisBase runs the
membrane-curvature analysis in every frame. We populate the arrays previously
defined in `_prepare()`. In `_conclude`, we compute the average over frames for
the surface and curvature arrays.

The derived surface, and obtained arrays of mean and Gaussian curvature values
are stored in the `Results` attribute. This makes of [AnalysisBase] the most
fundamental part of the MembraneCurvature analysis. With the membrane-curvature
AnalysisBase, we can perform multiframe and average-over-frame analysis.

We added also added coordinate wrapping to our Analysis base [#48], which
enables users to run MembraneCurvature with all atoms in the primary unit cell.
This is particularly useful when we want to calculate curvature in membrane-only
from a raw trajectory.

## Documentation ([#57], [#62], [#64], [#69])
One of the strongest motivations to contribute with a [MDAnalysis] 
curvature tool was to provide a well-documented package 
to analyze membrane curvature from MD simulations.

The membrane curvature tool includes solid documentation that can be found in the
following pages:

- [API documentation](https://membrane-curvature.readthedocs.io/en/latest/api.html) ([#57])
- [Algorithm], [Usage] and [Visualization] pages ([#62])
- [Tutorials] ([#64], [#69])

We included two different tutorials: One where we use Membrane Curvature to
derive surfaces and calculate curvature of a membrane-only system ([#64]), and
another for a membrane-protein system([#69]).


# How can I use it?
Membrane-curvature uses [MDAnalysis] under the hood. We can install Membrane-curvature via `pip`:

```
pip install membrane-curvature
```

MembraneCurvature was designed to be user friendly. No counterintuitive commands,
and no long lines of code. With MembraneCurvature you can calculate curvature in less
than 10 lines of code! The snippet below illustrates how easy it gets to extract mean
and Gaussian curvature from our MD simulations!


```
import MDAnalysis as mda
from membrane_curvature.base import MembraneCurvature
from membrane_curvature.tests.datafiles import MEMB_GRO, MEMB_XTC

u = mda.Universe(MEMB_GRO, MEMB_XTC)

curvature_upper_leaflet = MembraneCurvature(universe,
                                            select="resid 103-1023 and name PO4",
                                            n_x_bins=12,
                                            n_y_bins=12).run()


mean_upper_leaflet = curvature_upper_leaflet.results.average_mean
mean_lower_leaflet = curvature_lower_leaflet.results.average_mean

gaussian_upper_leaflet = curvature_upper_leaflet.results.average_gaussian
gaussian_lower_leaflet = curvature_lower_leaflet.results.average_gaussian
```

To visualize the results obtained in the analysis shown above, we can use
[contours] from matplotlib to get cool plots like the following:

```
from scipy import ndimage

curvatures = 
leaflets = ['Lower', 'Upper']
curvatures = [mean_lower_leaflet, mean_upper_leaflet]
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(5,3))
for ax, mc, lf in zip((ax1, ax2), curvatures, leaflets):
    ax.contourf(ndimage.gaussian_filter(mc, sigma=1, order=0, mode='reflect'),
                cmap='bwr',
                levels=30)
    ax.set_aspect('equal')
    ax.set_title('{} Leaflet'.format(lf))
    ax.axis('off')
```


![contours](../../assets/images/my_contours.png)


# In progress
Currently, we are working on implementing an interpolation as an option for the
user [#52].

In some situations, when selecting a very high number of bins in the grid, or
when having regions of the grid with low sampling. For example, think of a
membrane-protein system, where the bins occupied by the protein won't be
populated by lipids, and therefore, will have a region of undefined values in
the grid. Such undefined values spread in the array during the calculation of
curvature, which may result in meaningless results.

By adding an optional interpolation, we will be able to patch up undefined values 
in bins inside the embedded element (i.e. protein). With this improvement, calculation
of membrane curvature won't be hamstrung by the presence of undefined values in the grid.


# What's next?
There is always room for improvement, and MembraneCurvature is not an exception.

![vesicles](../../assets/images/vesicles.png)

One of the main limitations of the current version of MembraneCurvature is the
inability to calculate curvature in systems like vesicles, capsids, or micelles.
Definitely, this would be a nice improvement for a future release of
MembraneCurvature! We acknowledge that scientific research would benefit of a
tool to calculate membrane curvature in these types of systems, so we are
considering possible approaches!

# Conclusions
Participating in GSoC with MDAnalysis has been a unique experience. I had the
opportunity to learn best practices in software development mentored by a group
of incredibly talented people:
[@lilywang](https://github.com/lilyminium),
[@IAlibay](https://github.com/IAlibay), and
[@orbeckst](https://github.com/orbeckst). I also would like to thank
[@richardjgowers](https://github.com/richardjgowers) and 
[@tylerreddy](https://github.com/tylerjereddy) from the MDA community, who participated in our
discussions and provided valuable insights. 

Thanks for all your valuable lessons. 

MembraneCurvature has launched! 🚀

---
[#9]: https://github.com/MDAnalysis/membrane-curvature/pull/9
[#34]: https://github.com/MDAnalysis/membrane-curvature/pull/34
[#40]: https://github.com/MDAnalysis/membrane-curvature/pull/40
[#48]: https://github.com/MDAnalysis/membrane-curvature/pull/48
[#43]: https://github.com/MDAnalysis/membrane-curvature/pull/43
[#44]: https://github.com/MDAnalysis/membrane-curvature/pull/40
[#52]: https://github.com/MDAnalysis/membrane-curvature/pull/52
[#57]: https://github.com/MDAnalysis/membrane-curvature/pull/57
[#62]: https://github.com/MDAnalysis/membrane-curvature/pull/62
[#64]: https://github.com/MDAnalysis/membrane-curvature/pull/64
[#66]: https://github.com/MDAnalysis/membrane-curvature/pull/66
[#69]: https://github.com/MDAnalysis/membrane-curvature/pull/69

[MDAnalysis]: https://github.com/MDAnalysis
[Algorithm]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Algorithm.html
[Usage]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Usage.html
[Visualization]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Visualization.html
[Tutorials]: https://membrane-curvature.readthedocs.io/en/latest/source/pages/Tutorials.html
[AnalysisBase]: https://docs.mdanalysis.org/2.0.0-dev0/documentation_pages/analysis/base.html?highlight=base#analysis-building-blocks-mdanalysis-analysis-base