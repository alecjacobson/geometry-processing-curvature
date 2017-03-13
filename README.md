# Geometry Processing – Curvature 

> **To get started:** Fork this repository then issue
> 
>     git clone --recursive http://github.com/[username]/geometry-processing-curvature.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` by running
on a given mesh:

    ./curvature [path to mesh.obj]

## Background

In this assignment we explore discrete curvature quantities computed on a
surface. These quantities give us local information about a shape. Beyond
inspecting the surface (the extent of this assignment), these quantities become
the building blocks to:

 - define energies to minimize during smoothing/deformation,
 - identify salient points and curves on the shape, and 
 - provide initial conditions/constraints for remeshing.

The fundamental difference between a segment on the real line and a curve is
the introduction of 
[curvature](https://en.wikipedia.org/wiki/Curvature). This is quite natural and
intuitive. When we draw a 1D object in the plane or in space we have the
freedom to let that object bend. We quantify this "bending" locally as
curvature.

[Curvature](https://en.wikipedia.org/wiki/Curvature#Surfaces) is also the
fundamental difference between a chunk (i.e., subregion) of the [Euclidean
Plane](https://en.wikipedia.org/wiki/Plane_(geometry)) and a
[surface](https://en.wikipedia.org/wiki/Surface_(mathematics)) that has been
[immersed](https://en.wikipedia.org/wiki/Immersion_(mathematics)) in $\R^3$ (or
elsewhere). Unlike curves, surfaces can bend in each direction at any point.

## Tasks

- Gaussian curvature

- mean curvature 
