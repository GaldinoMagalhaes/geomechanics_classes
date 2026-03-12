# Rock Mechanics – Support Materials

## Initial Considerations

The materials shared here have been developed to support the undergraduate course Rock Mechanics at Polytechnique Montréal. The notebooks released so far focus on the topic of rock anisotropy and include:

- **Part I — Deformation of Anisotropic Rocks**
- **Part II — Anisotropic Failure Criteria**
- **Part III — Structure Tensors**

The materials are authored by the Geomechanics & Rock Mass Modeling research group and can be used for educational purposes, provided they are properly referenced. For suggestions, corrections, or general contact, please use the email [gabrielgaldinodm@gmail.com](mailto:gabrielgaldinodm@gmail.com). 

All shared content is referenced, and we highly recommend that it always be used in conjunction with the literature, as the primary focus of the material is conceptual review combined with programming exercises.

## Accessing the Notebooks

The Julia language was chosen for its ease in handling matrix data, high processing speed, and the interactive capabilities provided by the Pluto library. The notebooks can be accessed either statically via HTML files or interactively using the compiler. To acesses the interative notebooks, use the following steps:

1. Download the Julia compiler from the [Julia website](https://julialang.org/downloads/), or using the [MicrosoftStore](https://apps.microsoft.com/detail/9NJNWW8PVKMN?hl=en-us&gl=CA&ocid=pdpshare), for Windows users.

2. Download the [Notebooks](https://github.com/GaldinoMagalhaes/geomechanics_classes/tree/cfe29d918778ca5f234f176bb27b611cf6c561ef/notebooks) from the repository.

3. Install and Open Pluto, following:

A) Open the **Julia REPL**.  
B) Press `]` to enter the package manager (`(@v1.x) pkg>`).  
C) Install Pluto by typing:
```julia
add Pluto
```
D) Run Pluto by typing:
```julia
import Pluto
Pluto.run()
```
E) Input the notebooks path and click open.
