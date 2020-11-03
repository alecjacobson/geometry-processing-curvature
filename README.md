# Geometry Processing - Curvature 

> **To get started:** Clone this repository then issue
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

![Once built you may toggle between showing Gaussian curvature, mean/min/max
curvature and displaying principal curvature directions.](images/elephant-curvatures.jpg)

## Background

In this assignment we explore discrete curvature quantities computed on a
surface. These quantities give us local information about a shape. Beyond
inspecting the surface (the extent of this assignment), these quantities become
the building blocks to:

 - define energies to minimize during smoothing/deformation,
 - identify salient points and curves on the shape, and 
 - provide initial conditions/constraints for _remeshing_.

The fundamental difference between a segment on the real line and a curve is
the introduction of [curvature](https://en.wikipedia.org/wiki/Curvature). This
is quite natural and intuitive. When we draw a 1D object in the plane or in
space we have the freedom to let that object bend. We quantify this "bending"
locally as curvature.

[Curvature](https://en.wikipedia.org/wiki/Curvature#Surfaces) is also the
fundamental difference between a chunk (i.e., subregion) of the [Euclidean
Plane](https://en.wikipedia.org/wiki/Plane_(geometry)) and a
[surface](https://en.wikipedia.org/wiki/Surface_(mathematics)) that has been
[immersed](https://en.wikipedia.org/wiki/Immersion_(mathematics)) in <img alt="$\mathbb{R}^3$" src="./tex/d03c1e146df015e061405cc425738d83.svg" align="middle" width="18.424726649999986pt" height="26.76175259999998pt"/> (or
elsewhere). Unlike curves, surfaces can bend in each direction at any point.

We start our discussion assuming a smooth surface <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/>. We would like to
categorize points on the surface <img alt="$\mathbf{p} \in  \mathbf{S}$" src="./tex/5ba40130b56a36da7314110b06d636fe.svg" align="middle" width="41.09559299999999pt" height="22.55708729999998pt"/> in terms of how the surface bends or
curves locally. 

### Curvature of planar curves

Let us briefly recall how
[curvature](https://en.wikipedia.org/wiki/Curvature#Precise_definition) is
defined for a [planar curve](https://en.wikipedia.org/wiki/Plane_curve)
<img alt="${\gamma}:[0,1] \rightarrow  \mathbb{R}^{2}$" src="./tex/dd5ee1cd6ae3c66ca44f96af98021ce1.svg" align="middle" width="99.99435929999997pt" height="26.76175259999998pt"/>.

There are multiple equivalent definitions.

#### Osculating circle

We can define the [tangent](https://en.wikipedia.org/wiki/Tangent) direction at
a point <img alt="$\mathbf{p} = {\gamma}(s)$" src="./tex/a3aa09be6bf777e23fb57a2f3cb00887.svg" align="middle" width="62.334632249999984pt" height="24.65753399999998pt"/> as the limit of the
[secant](https://en.wikipedia.org/wiki/Secant_line) formed between <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> and
another point on the curve <img alt="$\mathbf{q} = {\gamma}(t)$" src="./tex/bb6c9ad7258cdd805398547a25e53aac.svg" align="middle" width="60.040138949999985pt" height="24.65753399999998pt"/> as <img alt="$\mathbf{q}$" src="./tex/e73485aa867794d51ccd8725055d03a3.svg" align="middle" width="9.97711604999999pt" height="14.611878600000017pt"/> approaches <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/>:

![](images/tangent-as-limit-of-secant.gif)

<p align="center"><img alt="$$&#10;\mathbf{t}(s) = \lim_{\mathbf{q}\rightarrow \mathbf{p}} \frac{\mathbf{q}-\mathbf{p}}{\|\mathbf{q}-\mathbf{p}\|} = &#10;\lim_{t\rightarrow s} \frac{{\gamma}(t)-{\gamma}(s)}{\|{\gamma}(t)-{\gamma}(s)\|} = \frac{{\gamma}'(s)}{\|{\gamma}'(s)\|}.&#10;$$" src="./tex/c001a1e15253b24d9c5c87222903ca82.svg" align="middle" width="371.54143949999997pt" height="38.864210549999996pt"/></p>

It always possible, and often convenient, to assume without loss of generality
that <img alt="$s$" src="./tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg" align="middle" width="7.7054801999999905pt" height="14.15524440000002pt"/> is an [arc length
parameterization](https://en.wikipedia.org/wiki/Arc_length) of the curve <img alt="${\gamma}$" src="./tex/193089f7a231633473714830d2edc62a.svg" align="middle" width="9.423880949999988pt" height="14.15524440000002pt"/> so
that <img alt="$\|{\gamma}'\| = 1$" src="./tex/41298c93d4b1388e9d8cf7017938b420.svg" align="middle" width="60.61099604999999pt" height="24.7161288pt"/> and therefor the unit tangent vector is simply <img alt="$\mathbf{t}(s) =&#10;{\gamma}'(s)$" src="./tex/5f8e197a7e88cc7cd47d4ed2a8f3ebd4.svg" align="middle" width="84.28675199999999pt" height="24.7161288pt"/>.

In an analogous fashion, we can consider the limit of the
[circumcircle](https://en.wikipedia.org/wiki/Circumscribed_circle)
<img alt="$C(\mathbf{q}_{1},\mathbf{p},\mathbf{q}_{2})$" src="./tex/4f92d02acce93762cc7226907c43394f.svg" align="middle" width="85.52721704999998pt" height="24.65753399999998pt"/> that passes
through <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> and points <img alt="$\mathbf{q}_{1}$" src="./tex/523c6393dc3eeaf6820cb573577cc8a0.svg" align="middle" width="16.529662049999992pt" height="14.611878600000017pt"/> and <img alt="$\mathbf{q}_{2}$" src="./tex/2643a8ebedff6c74098a107a85570ec4.svg" align="middle" width="16.529662049999992pt" height="14.611878600000017pt"/> before and after it on the curve:

<p align="center"><img alt="$$&#10;C(\mathbf{p}) = \lim_{\mathbf{q}_{1},\mathbf{q}_{2}\rightarrow \mathbf{p}} C(\mathbf{q}_{1},\mathbf{p},\mathbf{q}_{2}).&#10;$$" src="./tex/b368b0d08615562f235f543e4745cbbe.svg" align="middle" width="204.78055949999998pt" height="24.4292268pt"/></p>


![](images/osculating-circle.gif)

This limit circle is called the [osculating
circle](https://en.wikipedia.org/wiki/Osculating_circle) at the point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> on
the curve <img alt="${\gamma}$" src="./tex/193089f7a231633473714830d2edc62a.svg" align="middle" width="9.423880949999988pt" height="14.15524440000002pt"/>. By construction the tangent of the curve and the circle match at
<img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/>: they're both <img alt="${\gamma}'$" src="./tex/57c545ba93b947836e6461ddd4a3982f.svg" align="middle" width="13.213823699999992pt" height="24.7161288pt"/>. The
[radius](https://en.wikipedia.org/wiki/Radius_of_curvature) <img alt="$R(\mathbf{p})$" src="./tex/df2505b4579144078667b44d00513974.svg" align="middle" width="35.89613114999999pt" height="24.65753399999998pt"/> of the
osculating circle <img alt="$C(\mathbf{p})$" src="./tex/cf851d992ef882497c7bb5571aeb402a.svg" align="middle" width="36.21229754999999pt" height="24.65753399999998pt"/> at the the point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> is proportional to how straight
the curve is locally: as the curve becomes more and more straight then the
radius tends toward infinity. This implies that the radius is inversely
proportional to the "curvy-ness" of the curve. Hence, the inverse of the radius
is dubbed the curvature:

<p align="center"><img alt="$$&#10;{\kappa}(\mathbf{p}) = \frac{1}{R(\mathbf{p})}.&#10;$$" src="./tex/2de026fcd66c43e358680fd2f8a05401.svg" align="middle" width="99.08394704999999pt" height="37.099754999999995pt"/></p>

The radius is a non-negative measure of length with units meters, so the
curvature <img alt="${\kappa}$" src="./tex/2af59c6d7260c35c9ddd118c6d40f5d9.svg" align="middle" width="9.47111549999999pt" height="14.15524440000002pt"/> is an non-negative scalar with units 1/meters. The radius of the
osculating circle can also be written as a limit of the [circumcircle
radius](https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_from_cross-_and_dot-products):

<p align="center"><img alt="$$&#10;R(\mathbf{p}) = \lim_{\mathbf{q}_{1},\mathbf{q}_{2}\rightarrow \mathbf{p}} &#10;  \frac{\| \mathbf{q}_{1}-\mathbf{p}\|  \| \mathbf{p}-\mathbf{q}_{2}\|  \| \mathbf{q}_{2}-\mathbf{q}_{1}\| }&#10;  {2\left| (\mathbf{q}_{1}-\mathbf{p}) \quad (\mathbf{p}-\mathbf{q}_{2})\right|}.&#10;$$" src="./tex/06eac18ac26c6374d00e5b44f56aa39b.svg" align="middle" width="322.8819033pt" height="38.83491479999999pt"/></p>


#### Signed curvature

Plugging in our arc-length parameterization this reveals that the curvature
(inverse of radius) is equal to the magnitude of change in the tangent or
equivalently the magnitude of second derivative of the curve:

<p align="center"><img alt="$$&#10;{\kappa}(s) = \lim_{t\rightarrow s} \left\| \frac{{\gamma}'(t)-{\gamma}'(s)}{t-s} \right| = \| {\gamma}''(s)\| .&#10;$$" src="./tex/b599230eafe7e17b9040c86efd6ef708.svg" align="middle" width="266.06554919999996pt" height="39.452455349999994pt"/></p>


Because we chose the arc-length parameterization, the only change to the
tangent vector <img alt="${\gamma}'$" src="./tex/57c545ba93b947836e6461ddd4a3982f.svg" align="middle" width="13.213823699999992pt" height="24.7161288pt"/> is a change in _direction_ (as opposed to magnitude, since
<img alt="$\| {\gamma}'\|  := 1$" src="./tex/fc7bc861574810d4b6da99ff1c77f1a2.svg" align="middle" width="65.17721924999998pt" height="24.7161288pt"/>). This means that the change--as a vector itself--is
_orthogonal_ to the tangent. In other words, the change in tangent <img alt="${\gamma}''$" src="./tex/387e1fec91a9dec19df1c99c8ecd15bf.svg" align="middle" width="17.003784599999992pt" height="24.7161288pt"/> points
along the normal direction <img alt="$\widehat{\mathbf{n}}$" src="./tex/65864fff030088096bc2384e9712668f.svg" align="middle" width="10.502226899999991pt" height="24.200985600000003pt"/>:

<p align="center"><img alt="$$&#10;{\gamma}'' \cdot  {\gamma}' = 0 \quad \rightarrow  \quad {\gamma}'' \cdot  \widehat{\mathbf{n}} = \pm  {\kappa} \widehat{\mathbf{n}}.&#10;\label{equ:curvature-normal}&#10;$$" src="./tex/fdb8daaccc56847b848b7a8ce1c2f63f.svg" align="middle" width="231.76014674999996pt" height="16.3763325pt"/></p>


If we define an orientation to our curve  then we can endow the curvature with
a [sign](https://en.wikipedia.org/wiki/Sign_(mathematics)) based on whether the
center of the osculating circle lies on the [left or right
side](https://en.wikipedia.org/wiki/Right-hand_rule) of the curve. As already
established, the tangent of the osculating circle and the curve agree, so the
vector pointing toward the circle's center must be
[perpendicular](https://en.wikipedia.org/wiki/Perpendicular) to the tangent:
i.e., in either the positive or negative
[normal](https://en.wikipedia.org/wiki/Normal_(geometry)) directions. 

If the orientation agrees with increasing the arc-length parameter <img alt="$s$" src="./tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg" align="middle" width="7.7054801999999905pt" height="14.15524440000002pt"/>, then the sign can
be 
determined by comparing the second derivative vector <img alt="${\gamma}''$" src="./tex/387e1fec91a9dec19df1c99c8ecd15bf.svg" align="middle" width="17.003784599999992pt" height="24.7161288pt"/> to the unit normal
<img alt="$\widehat{\mathbf{n}} :=&#10;({\gamma}')^{\perp}$" src="./tex/76ec459449f816b5c029cc27353af134.svg" align="middle" width="74.08127264999999pt" height="27.91243950000002pt"/>. The [_**signed
curvature**_](https://en.wikipedia.org/wiki/Curvature#Signed_curvature) at a point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> is thus given by:

<p align="center"><img alt="\begin{align*}&#10;k(\mathbf{p}) &amp;= \text{sign}({\gamma}''(\mathbf{p})\cdot \widehat{\mathbf{n}}))\ {\kappa}(\mathbf{p})  \\&#10;      &amp;= {\gamma}''(\mathbf{p}) \cdot  \widehat{\mathbf{n}}.&#10;\end{align*}" src="./tex/03bdb9a620549fa6bdb44f036303d7e7.svg" align="middle" width="203.5864677pt" height="41.9471052pt"/></p>


#### Moving point analogy

This definition neatly conforms to our intuition of a curve as the trajectory
of a moving point. Imagine the curved formed by driving along a particular
trajectory <img alt="${\gamma}(t)$" src="./tex/9f34318b41b902d4f823fe84b7c99d54.svg" align="middle" width="28.14539309999999pt" height="24.65753399999998pt"/>, where we really interpret <img alt="$t$" src="./tex/4f4f4e395762a3af4575de74c019ebb5.svg" align="middle" width="5.936097749999991pt" height="20.221802699999984pt"/> as time.

While <img alt="${\gamma}'(t)$" src="./tex/ec3c193e6dced0aef0aa4a5f211a09c2.svg" align="middle" width="32.757266849999986pt" height="24.7161288pt"/> corresponds to your velocity vector and <img alt="$\| {\gamma}'(t)\| $" src="./tex/8f893124cc00bec18bb594971331a09a.svg" align="middle" width="49.195685549999986pt" height="24.7161288pt"/> corresponds to
your speed, the arc-length (re-)parameterization would correspond to having
your friend re-trace your path traveling at a perfectly uniform speed <img alt="$\| {\gamma}'(s)\| &#10;= 1$" src="./tex/b24cdd55db07f787cd207004056ec33a.svg" align="middle" width="81.10190879999999pt" height="24.7161288pt"/>, where your friends "time" <img alt="$s$" src="./tex/6f9bad7347b91ceebebd3ad7e6f6f2d1.svg" align="middle" width="7.7054801999999905pt" height="14.15524440000002pt"/> may be different from yours (it may take
longer or shorter depending if you drove fast or slow).

Curvature in the path corresponds to _turning_ and quite literally the amount
by which your friend needs to turn the steering wheel away from the "straight"
position: on a straight course, the steering wheel remains at zero-angle
position and the curvature is zero, on a circular course the steering wheel is
fixed at a constant angle in the left or right direction corresponding to
constant positive or negative curvature respectively.

Changing the steering wheel changes the _direction_ of the vehicle's velocity.
For your friend driving at constant speed, this is the _only_ change admissible
to the velocity, hence the curvature exactly corresponds to <img alt="${\gamma}''(s)$" src="./tex/e9600f89417cea0ec19d927bc0825861.svg" align="middle" width="38.31660689999999pt" height="24.7161288pt"/> and to the
steering wheel angle.

> If somebody wants to make a Sega [Out
> Run](https://en.wikipedia.org/wiki/Out_Run) inspired gif showing a steering
> wheel turning next to a little car tracing a curve, I'll be very impressed.

#### Turning number

The integrated signed curvature around a [closed
curve](https://en.wikipedia.org/wiki/Curve) must be an integer multiple of
<img alt="$2{\pi}$" src="./tex/a04c31f698b09c6dd0c7687b3514164e.svg" align="middle" width="18.179315549999988pt" height="21.18721440000001pt"/>:

<p align="center"><img alt="$$&#10;{\oint} k(s) \ ds = 2{\pi} {\tau},&#10;$$" src="./tex/439b2ad254f6069803a4745194ef4f14.svg" align="middle" width="120.5422845pt" height="36.53007435pt"/></p>

where <img alt="$\tau$" src="./tex/0fe1677705e987cac4f589ed600aa6b3.svg" align="middle" width="9.046852649999991pt" height="14.15524440000002pt"/> is an integer called the "turning number" of the curve. 

This is a bit surprising at first glance. However, in the _moving point
analogy_ a closed curve corresponds to a period trajectory (e.g., driving
around a race-track). When we've made it once around the track, our velocity
direction (e.g., the direction the vehicle is facing) must be pointing in the
original direction. That is, during the course, the car either have turned all the
way around once (<img alt="${\tau} = 1$" src="./tex/2e382d1d2f6b4ef71deca3309c5f80f0.svg" align="middle" width="39.18367364999999pt" height="21.18721440000001pt"/>) or turned as much clockwise and it did
counter-clockwise (e.g., on a figure 8 course: <img alt="${\tau}=0$" src="./tex/5caa88446701900327eb3af4dbbf5126.svg" align="middle" width="39.18367364999999pt" height="21.18721440000001pt"/>), or made multiple
loops, etc.

#### Discrete curvature

In the discrete world, if a curve is represented as a piecewise-linear chain of
segments, then it's natural to associate curvature with vertices: the segments
are flat and therefor contain no curvature.

A natural analog to the definition of curvature as 
the derivative of the tangent vector 
(i.e., <img alt="$k = \| {\gamma}''\|  = \| \mathbf{t}'\| $" src="./tex/9f3507bcfeb67b4a30abfddfc3aee597.svg" align="middle" width="115.57658969999997pt" height="24.7161288pt"/>) is to define _discrete curvature_ as the change in
tangent direction between discrete segments meeting at a vertex:

<p align="center"><img alt="$$&#10;k_i =  {\angle} (\mathbf{v}_i - (\mathbf{v}_{i-1}-\mathbf{v}_i)) \mathbf{v}_i\mathbf{v}_{i+1} = {\theta}_i,&#10;$$" src="./tex/6473c8206b40ed854546bb8298ec8dc6.svg" align="middle" width="263.96728874999997pt" height="16.438356pt"/></p>

that is, the signed [_exterior
angle_](https://en.wikipedia.org/wiki/Internal_and_external_angles) <img alt="${\theta}_i$" src="./tex/be9fa78ef5e17dc25159197f21d1e3cb.svg" align="middle" width="12.36779114999999pt" height="22.831056599999986pt"/> at
the vertex <img alt="$\mathbf{v}_i$" src="./tex/5474c8baa2bc6feabb0eac4237772aab.svg" align="middle" width="14.628015599999989pt" height="14.611878600000017pt"/>.

![](images/discrete-curvature.svg)

The turning number theorem for continuous curves finds an _immediate_ analog in
the discrete case. For a closed polygon the discrete signed angles must sum up
to a multiple of <img alt="$2{\pi}$" src="./tex/a04c31f698b09c6dd0c7687b3514164e.svg" align="middle" width="18.179315549999988pt" height="21.18721440000001pt"/> in order to close up:

<p align="center"><img alt="$$&#10;{\sum}_{i=1}^n k_i = 2{\pi} {\tau}.&#10;$$" src="./tex/bcaa4d69bdd3a84cf27ba108aa705001.svg" align="middle" width="112.68846104999999pt" height="28.013911200000003pt"/></p>


In this way, we _preserve the structure_ found in the continuous case in our
discrete analog. This structure preservation leads to an understanding of the
exterior angle as an approximation or discrete analog of the _locally
integrated_ curvature.

Alternatively, we could literally fit an circle to the discrete curve based on
local samples and approximate curvature as the inverse radius of the osculating
circle. This curvature measure (in general) will not obey the turning number
theorem, but (conducted properly) it will converge to the pointwise continuous
values under refinement (e.g., as segment length shrinks).

We will explore these two concepts for surfaces, too: discrete analogs that
preserve continuous structures and discretizations that approximate continuous
quantities in the limit.

### Curvature(s) on surfaces

A surface can be curved locally in multiple ways. Consider the difference
between a flat piece of paper, a spherical ping-pong ball and a saddle-shaped
[Pringles chip](https://en.wikipedia.org/wiki/Pringles). The Pringles chip is
the most interesting because it curves "outward" in one direction and "inward"
in another direction. In this section, we will learn to distinguish and
classify points on a surface based on how it curves in each direction.

#### Normal curvature

The simplest way to extend the curvature that we defined for planar curves to a
surface <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/> is to _slice_ the surface through a given point <img alt="$\mathbf{p}\in \mathbf{S}$" src="./tex/f6f3f7d110e410b0f3c012aca3c16b07.svg" align="middle" width="41.09559299999999pt" height="22.55708729999998pt"/> with a
[plane](https://en.wikipedia.org/wiki/Plane_(geometry)) <img alt="$\mathbf{P}$" src="./tex/384591906555413c452c93e493b2d4ec.svg" align="middle" width="12.92230829999999pt" height="22.55708729999998pt"/> that is parallel
to the [surface normal](https://en.wikipedia.org/wiki/Normal_(geometry))
<img alt="$\mathbf{n}(\mathbf{p})$" src="./tex/9ab0e9039036911be55bc8f0dc28ba8c.svg" align="middle" width="33.78988799999999pt" height="24.65753399999998pt"/>.

The (local) intersection of the surface <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/> and the plane <img alt="$\mathbf{P}$" src="./tex/384591906555413c452c93e493b2d4ec.svg" align="middle" width="12.92230829999999pt" height="22.55708729999998pt"/> will trace a
curve <img alt="${\gamma}$" src="./tex/193089f7a231633473714830d2edc62a.svg" align="middle" width="9.423880949999988pt" height="14.15524440000002pt"/>, upon which we can immediately use the planar curvature definition
above. 

![[Source](http://brickisland.net/cs177fa12/?p=214)](images/normal-curvature.svg)

There are infinitely many planes that pass through a given point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> and lie
parallel to a given normal vector <img alt="$\mathbf{n}(\mathbf{p})$" src="./tex/9ab0e9039036911be55bc8f0dc28ba8c.svg" align="middle" width="33.78988799999999pt" height="24.65753399999998pt"/>: the plane can rotate around the
normal <img alt="$\mathbf{n}(\mathbf{p})$" src="./tex/9ab0e9039036911be55bc8f0dc28ba8c.svg" align="middle" width="33.78988799999999pt" height="24.65753399999998pt"/> by any angle <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/>. For each choice of <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/>, the plane will define
an intersecting curve <img alt="${\gamma}_{\varphi}$" src="./tex/bbf36b32d259b8e3d56e0ecdee517461.svg" align="middle" width="17.15681219999999pt" height="14.15524440000002pt"/> and thus for every angle <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/> there will be a
[_normal curvature_](https://en.wikipedia.org/wiki/Curvature#Normal_sections):

<p align="center"><img alt="$$&#10;k_\mathbf{n}({\varphi},\mathbf{p}) = {\gamma}''_\mathbf{n}(\mathbf{p}).&#10;$$" src="./tex/7b2e2457fbcaddc566fc9b9dbca30c53.svg" align="middle" width="126.58880849999998pt" height="17.2895712pt"/></p>


#### Mean curvature

Normal curvature requires choosing an angle, so it doesn't satiate our
desire to reduce the "curvy-ness" to a single number for any point on the
surface. A simple way to reduce this space of normal curvatures is to, well,
average all possible normal curvatures. This defines the [mean
curvature](https://en.wikipedia.org/wiki/Mean_curvature):

<p align="center"><img alt="$$&#10;H(\mathbf{p}) = \frac{1}{2{\pi}}\int _0^{2{\pi}} k_\mathbf{n}({\varphi},\mathbf{p}) \ d{\varphi}.&#10;$$" src="./tex/33bc70f308d152540117375865d2fbe9.svg" align="middle" width="208.06852604999997pt" height="40.70359755pt"/></p>


#### Maximum and minimum curvature 

Another obvious way to reduce the space of normal curvatures to a single number
is to consider the maximum or minimum normal curvature over all choices of <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/>:

<p align="center"><img alt="\begin{align*}&#10;k_{1}(\mathbf{p}) &amp;= \max_{\varphi} \ k_\mathbf{n}({\varphi},\mathbf{p}) \\&#10;k_{2}(\mathbf{p}) &amp;= \mathop{\text{min}}_{\varphi} \ k_\mathbf{n}({\varphi},\mathbf{p}).&#10;\end{align*}" src="./tex/5d897003d90cebd14279c85c6bc39422.svg" align="middle" width="160.31188799999998pt" height="58.55700015pt"/></p>


Collectively, these are referred to as the [principal
curvatures](https://en.wikipedia.org/wiki/Principal_curvature) and
correspondingly the angles that maximize and minimize curvature are referred to
as the principal curvature directions:

<p align="center"><img alt="\begin{align*}&#10;{\varphi}_{1}(\mathbf{p}) &amp;= \mathop{\text{argmax}}_{\varphi} \ k_\mathbf{n}({\varphi},\mathbf{p}) \\&#10;{\varphi}_{2}(\mathbf{p}) &amp;= \mathop{\text{argmin}}_{\varphi} \ k_\mathbf{n}({\varphi},\mathbf{p}).&#10;\end{align*}" src="./tex/cd6fb93c75f85ab6cdc57cca14f9b0b7.svg" align="middle" width="185.38421879999999pt" height="64.9496892pt"/></p>


[Euler's
theorem](https://en.wikipedia.org/wiki/Euler%27s_theorem_(differential_geometry))
states that the normal curvature is a quite simple function of <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/> and the
principal curvatures:

<p align="center"><img alt="$$&#10;k_\mathbf{n}({\varphi},\mathbf{p}) = k_{1} \cos^2 {\varphi} + k_{2} \sin^2 {\varphi},&#10;$$" src="./tex/d1e44843b9156c97c125ab06c1cfac19.svg" align="middle" width="226.83757964999998pt" height="18.439728449999997pt"/></p>

([proof](https://math.stackexchange.com/a/1783316/35376)).

There are two immediate and important consequences:

 1. the principal curvature directions (<img alt="${\varphi}_{1}$" src="./tex/67ef67c97ad06a0f412383f1a0f35dd8.svg" align="middle" width="17.30598044999999pt" height="14.15524440000002pt"/> and <img alt="${\varphi}_{2}$" src="./tex/5fa5023351dd2803361b6b11163f1307.svg" align="middle" width="17.30598044999999pt" height="14.15524440000002pt"/>) are orthogonal, and 
 2. the mean curvature reduces to the average of principal curvatures:

<p align="center"><img alt="$$&#10;H = \frac12 (k_{1} + k_{2}).&#10;$$" src="./tex/bcf370af6bd62035c375d12e082f4248.svg" align="middle" width="118.38951794999998pt" height="32.990165999999995pt"/></p>


> For more theory and a proof of Euler's theorem, I recommend "Elementary
> Differential Geometry" by Barret O'Neill, Chapter 5.2.

#### Gaussian curvature

Maximum, minimum and mean curvature reduce curvature to a single number, but
still cannot (alone) distinguish between points lying on a round ping-pong
ball, a flat sheet of paper, the cylindrical Pringles can and a saddle-shaped
Pringles chip.

The neck of this cartoon elephant--like a Pringles chip--bends inward in one
direction (positive <img alt="$k_{1} &gt; 0$" src="./tex/ed6d6e66630f757e45b8b30927d3388c.svg" align="middle" width="46.069176449999986pt" height="22.831056599999986pt"/>) and outward in the other 
direction (negative <img alt="$k_{2} &lt; 0$" src="./tex/bf5d96e517812c3a5e6944e787aa4d5f.svg" align="middle" width="46.069176449999986pt" height="22.831056599999986pt"/>).

![Maximum <img alt="$k_{1}$" src="./tex/df2a706e439b4632680202a540f27fd0.svg" align="middle" width="15.11042279999999pt" height="22.831056599999986pt"/>, minimum <img alt="$k_{2}$" src="./tex/bce145d2ce5a01a17e36278910bfa8a4.svg" align="middle" width="15.11042279999999pt" height="22.831056599999986pt"/>,and Gaussian curvature <img alt="$K = k_{1}k_{2}$" src="./tex/f6fbf53c0380570ad3973fcc160e79c1.svg" align="middle" width="68.09738924999999pt" height="22.831056599999986pt"/>](images/cartoon-elephant-principal-and-gaussian-curvature.jpg)

The _product_ of the principal curvatures maintains the disagreement in sign
that categories this saddle-like behavior. This product is called [Gaussian
curvature](https://en.wikipedia.org/wiki/Gaussian_curvature):

<p align="center"><img alt="$$&#10;K = k_1 k_2.&#10;$$" src="./tex/1405680900c15675bf1512760b378718.svg" align="middle" width="73.48552695pt" height="13.881256950000001pt"/></p>


#### Relationship to surface area

Both mean and Gaussian curvature have meaningful relationships to surface
area.

##### Mean Curvature as area gradient

Let us consider a seemingly unrelated yet familiar problem. Suppose we would
like to minimize the surface area of a given shape <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/>. 

If we consider a small patch <img alt="$A \subset  \mathbf{S}$" src="./tex/99d36d49d13126eff89c1eccc623a5cf.svg" align="middle" width="44.74865669999999pt" height="22.55708729999998pt"/> centered around some point <img alt="$\mathbf{p} \in \mathbf{S}$" src="./tex/9f0b9d6cf50d1afd6a6923bfbf3b06a2.svg" align="middle" width="41.09559299999999pt" height="22.55708729999998pt"/>, then
we can imagine forces at the boundary of the patch <img alt="$\partial A$" src="./tex/539792c0748b098f331fc66c50e87231.svg" align="middle" width="21.96920384999999pt" height="22.831056599999986pt"/> that try to shrink the
area of the patch. These "surface tension" forces point inward against the
normal direction of the patch boundary <img alt="${\eta}$" src="./tex/7f2884b982daf1ca0cac1cf81bfeceec.svg" align="middle" width="8.751954749999989pt" height="14.15524440000002pt"/>. Consider the integral effect of
these forces and apply [divergence
theorem](https://en.wikipedia.org/wiki/Divergence_theorem) to move the integral
to the interior of the patch.

<p align="center"><img alt="$$&#10;{\oint}_{\partial A} -{\eta} \;ds = {\oint}_{\partial A} \mathbf{n} \times  d\mathbf{x} = \int _A \Delta \mathbf{x} \;dA,&#10;$$" src="./tex/6f3f1f19dd82f297215d5f3c3ddc0423.svg" align="middle" width="288.65304764999996pt" height="37.3519608pt"/></p>

where <img alt="$\mathbf{n}$" src="./tex/b56595d2a30a0af329086562ca12d521.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> is the surface normal, <img alt="$d\mathbf{x}$" src="./tex/f65261255ac2de833443ba78d717edb3.svg" align="middle" width="18.533080499999986pt" height="22.831056599999986pt"/> is the infinitesimal change in position
around the boundary of the patch, and <img alt="$\Delta \mathbf{x} \in  \mathbb{R}^{3}$" src="./tex/3af81c748738bde7407bd6771090767f.svg" align="middle" width="62.191651499999985pt" height="26.76175259999998pt"/> computes the
[Laplacian](https://en.wikipedia.org/wiki/Laplace_operator) of the position
coordinate functions.

In the limit as <img alt="$A$" src="./tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg" align="middle" width="12.32879834999999pt" height="22.465723500000017pt"/> shrinks around any point on the surface,
this relationship gives us the intuition that _moving_ the surface in the
direction <img alt="$\Delta \mathbf{x}$" src="./tex/52cdbba1eb4fddc47169867a1330fa0d.svg" align="middle" width="23.67578729999999pt" height="22.465723500000017pt"/> will decrease the surface area. 

The Laplacian <img alt="$\Delta f$" src="./tex/5d287dab9fa370f57a1a7ca52c46f728.svg" align="middle" width="23.516088749999987pt" height="22.831056599999986pt"/> of a function <img alt="$f$" src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/> on the surface does not depend on the
choice of parameterization. It is defined as the divergence of the gradient of
the function or equivalently the trace of the Hessian:

<p align="center"><img alt="$$&#10;\Delta f = {\nabla}\cdot  {\nabla}f = \tr{ &#10;\left[ &#10;\begin{array}{cc}&#10;\frac{\partial ^{2}f}{\partial u^{2}} &amp; \frac{\partial ^{2}f}{\partial u\partial v} \\&#10;\frac{\partial ^{2}f}{\partial v\partial u} &amp; \frac{\partial ^{2}f}{\partial v^{2}} &#10;\end{array}&#10;\right]&#10;}&#10;=&#10;\frac{\partial ^{2}f}{\partial u^{2}} + \frac{\partial ^{2}f}{\partial v^{2}}.&#10;$$" src="./tex/880dd562c4ea52694c6ec8007b6344a9.svg" align="middle" width="344.94506805pt" height="49.315569599999996pt"/></p>


If we generously choose <img alt="$u$" src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/> and <img alt="$v$" src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/> to vary in the principal directions <img alt="${\varphi}_{1}$" src="./tex/67ef67c97ad06a0f412383f1a0f35dd8.svg" align="middle" width="17.30598044999999pt" height="14.15524440000002pt"/>
and <img alt="${\varphi}_{2}$" src="./tex/5fa5023351dd2803361b6b11163f1307.svg" align="middle" width="17.30598044999999pt" height="14.15524440000002pt"/> above. In this case, the Laplacian <img alt="$\Delta \mathbf{x}$" src="./tex/52cdbba1eb4fddc47169867a1330fa0d.svg" align="middle" width="23.67578729999999pt" height="22.465723500000017pt"/> of the position function
reduces to the sum of principal curvatures times the normal (recall Equation
<img alt="$(\ref{equ:curvature-normal})$" src="./tex/947d5dbf6ad411cb13fe9fdb61a20183.svg" align="middle" width="30.63921959999999pt" height="24.65753399999998pt"/>):

<p align="center"><img alt="\begin{align*}&#10;\Delta \mathbf{x} &amp;= \frac{\partial ^{2}\mathbf{x}}{\partial u^{2}} + \frac{\partial ^{2}\mathbf{x}}{\partial v^{2}} \\&#10;    &amp;= k_{1} \mathbf{n} + k_{2} \mathbf{n} \\&#10;    &amp;= 2H\mathbf{n},&#10;\end{align*}" src="./tex/4bf6b8f7cfd3198e949aad1bc052fe3e.svg" align="middle" width="125.58636255pt" height="84.0148452pt"/></p>

where <img alt="$H\mathbf{n} \in  \mathbb{R}^{3}$" src="./tex/efcabc28b96b909eb308809defbb1867.svg" align="middle" width="64.01805959999999pt" height="26.76175259999998pt"/> is called the _**mean curvature normal**_ vector. We have
shown that the mean curvature normal is equal half the Laplacian of the
embedding function.

<!--
This turns out to be surprisingly complicated to derive without "going of the
deep end" into differential forms etc.

The trivial definition of surface area as an integral over <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/>
becomes an integral over <img alt="${\Omega}$" src="./tex/9f531c9f3f1ebeef802ced46eabb0336.svg" align="middle" width="11.87217899999999pt" height="22.465723500000017pt"/> of the determinant of the
[Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of
<img alt="$\mathbf{x}$" src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg" align="middle" width="9.97711604999999pt" height="14.611878600000017pt"/> after applying [substitution
rule](https://en.wikipedia.org/wiki/Integration_by_substitution):

<p align="center"><img alt="$$&#10;\int _\mathbf{S} 1\ dA = \int _{\Omega} \left| {\nabla}\mathbf{x} \right| d{\Omega}&#10;$$" src="./tex/2aa4d55c21ff67ec8ab59fff816d9abb.svg" align="middle" width="155.5135659pt" height="37.3519608pt"/></p>


-->

##### Gaussian Curvature as area distortion

As the product of principal curvatures, Gaussian curvature <img alt="$K = k_{1}k_{2}$" src="./tex/f6fbf53c0380570ad3973fcc160e79c1.svg" align="middle" width="68.09738924999999pt" height="22.831056599999986pt"/> measures zero
anytime one (or both) of the principal curvatures are zero. Intuitively, this
happens only for surfaces that curve or bend in one direction. Imagine rolling
up a sheet of paper. Surfaces with zero Gaussian curvature <img alt="$K = 0$" src="./tex/a353f3e8572d45386453eaa4b27f0bc2.svg" align="middle" width="45.273840149999984pt" height="22.465723500000017pt"/> are called
_developable surfaces_ because the can be flattened (developed) on to the flat
plane (just as you might unroll the piece of paper) _without_ stretching or
shearing. As a corollary, surfaces with non-zero Gaussian curvature _cannot_ be
flattened to the plane without stretching some part.

Locally, Gaussian curvature measures how far from developable the surface is:
how much would the local area need to stretch to become flat.

First, we introduce the [Gauss map](https://en.wikipedia.org/wiki/Gauss_map), a
continuous map <img alt="$N:\mathbf{S}\rightarrow S^{2}$" src="./tex/010cd709f53ffd8c3d7357e2eb093ffa.svg" align="middle" width="82.35114689999999pt" height="26.76175259999998pt"/> from every point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> on the surface <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/> to the unit
sphere <img alt="$S^{2}$" src="./tex/d95580c1170185b7f1bf99f798ef9fcc.svg" align="middle" width="17.57992994999999pt" height="26.76175259999998pt"/> so that <img alt="$N(\mathbf{p}) := \mathbf{n}(\mathbf{p})$" src="./tex/847e60608c79423b277b2e7bd9e6bf34.svg" align="middle" width="98.56137224999998pt" height="24.65753399999998pt"/>, the unit normal at <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/>.

Consider a small patch on a curved surface. Gaussian curvature <img alt="$K$" src="./tex/d6328eaebbcd5c358f426dbea4bdbf70.svg" align="middle" width="15.13700594999999pt" height="22.465723500000017pt"/> can
equivalently be defined as the limit of the ratio between the area
area _swept_ out by the unit normal on the Gauss map <img alt="$A_G$" src="./tex/85b403df61d71bf1dd88dbaf8579deb9.svg" align="middle" width="22.563292949999987pt" height="22.465723500000017pt"/> and 
the area of the surface patch <img alt="$A$" src="./tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg" align="middle" width="12.32879834999999pt" height="22.465723500000017pt"/>:

<p align="center"><img alt="$$&#10;K = \lim_{A\rightarrow 0} \frac{A_G}{A}.&#10;\label{equ:gaussian-curvature-area}&#10;$$" src="./tex/d2167c751819a227677e0ab1c5be9219.svg" align="middle" width="101.23440194999999pt" height="33.62942055pt"/></p>


Let's consider different types of regions:

 - flat: <img alt="$A_G=0$" src="./tex/c7fd05a94bba5f33ad41ca4789012009.svg" align="middle" width="53.52204824999999pt" height="22.465723500000017pt"/> because the Gauss map is a point,
 - cylindrical: <img alt="$A_G=0$" src="./tex/c7fd05a94bba5f33ad41ca4789012009.svg" align="middle" width="53.52204824999999pt" height="22.465723500000017pt"/> because the Gauss map is a curve,
 - spherical: <img alt="$A_G &gt; 0$" src="./tex/2d1cf747af9ff97e3754507d60c8c75a.svg" align="middle" width="53.52204824999999pt" height="22.465723500000017pt"/> because the Gauss map will maintain positive swept-area,
   and 
 - saddle-shaped: <img alt="$A_G &lt; 0$" src="./tex/b1263473b0f4fd328120721193ffeb07.svg" align="middle" width="53.52204824999999pt" height="22.465723500000017pt"/> because the area on the Gauss map will maintain
   _oppositely_ oriented area (i.e., from the spherical case).

![A patch on a plane and its corresponding patch on the Gauss map.](images/gauss-map-plane.gif)
![A patch on a cylinder and its corresponding patch on the Gauss map.](images/gauss-map-cylinder.gif)
![A patch on a sphere and its corresponding patch on the Gauss map.](images/gauss-map-sphere.gif)
![A patch on a saddle and its corresponding patch on the Gauss map.](images/gauss-map-saddle.gif)

Similar to the turning number theorem for curves, there exists an analogous
[theorem for surfaces](https://en.wikipedia.org/wiki/Gauss-Bonnet_theorem)
stating that the total Gaussian curvature must be an integer multiple of <img alt="$2{\pi}$" src="./tex/a04c31f698b09c6dd0c7687b3514164e.svg" align="middle" width="18.179315549999988pt" height="21.18721440000001pt"/>:

<p align="center"><img alt="$$&#10;\int _S K dA = 2{\pi} {\chi}(\mathbf{S}),&#10;\label{equ:gauss-bonnet}&#10;$$" src="./tex/17545b2bf4691f8949e50b818a2be4fb.svg" align="middle" width="135.65290035pt" height="37.3519608pt"/></p>

where <img alt="${\chi}(\mathbf{S})$" src="./tex/3d04c4e853af7437cfe553ef3343c970.svg" align="middle" width="33.57301529999999pt" height="24.65753399999998pt"/> is the [Euler
characteristic](https://en.wikipedia.org/wiki/Euler_characteristic) of the
surfaces <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/> (a topological _invariant_ of the surface revealing how many
[holes](https://en.wikipedia.org/wiki/Genus_(mathematics)) the surface has).

In stark contrast to mean curvature, this theorem tells us that we cannot add
Gaussian curvature to a surface without:

  1. removing an equal amount some place else, or
  2. changing the topology of the surface.

Since changing the topology of the surface would require a discontinuous
deformation, adding and removing Gaussian curvature must also balance out for
smooth deformations. This simultaneously explains why a cloth must have
wrinkles when draping over a table, and why a deflated basketball will not lie
flat on the ground.

#### Shape operator

There is yet another way to arrive at principal, mean and Gaussian curvatures.
Consider a point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> on a surface <img alt="$\mathbf{S}$" src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg" align="middle" width="10.502226899999991pt" height="22.55708729999998pt"/> with unit normal vector <img alt="$\mathbf{n}$" src="./tex/b56595d2a30a0af329086562ca12d521.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/>. If we
pick a unit tangent vector <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> (i.e., so that <img alt="$\mathbf{v} \cdot  \mathbf{n} = 0$" src="./tex/db0d91aa7e9eb5266fcc930eb1b02f4b.svg" align="middle" width="62.75073254999998pt" height="21.18721440000001pt"/>), then we can ask
how does the normal <img alt="$\mathbf{n}$" src="./tex/b56595d2a30a0af329086562ca12d521.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> change as we move in the direction of <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> along the
surface:

<p align="center"><img alt="$$&#10;S_\mathbf{p}(\mathbf{v}) := {\nabla} \mathbf{n} \cdot  \mathbf{v}&#10;$$" src="./tex/f066356ca1ab074b3dad6bb398a824b8.svg" align="middle" width="114.98821124999999pt" height="17.031940199999998pt"/></p>

we call <img alt="$S_\mathbf{p}$" src="./tex/2296c3cd641c58ee180596a79543cef2.svg" align="middle" width="18.34475609999999pt" height="22.465723500000017pt"/> the [_**shape
operator**_](https://en.wikipedia.org/wiki/Differential_geometry_of_surfaces#Shape_operator)
at the point <img alt="$\mathbf{p}$" src="./tex/980fcd4213d7b5d2ffcc82ec78c27ead.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/>. Just as how in Equation <img alt="$(\ref{equ:curvature-normal})$" src="./tex/947d5dbf6ad411cb13fe9fdb61a20183.svg" align="middle" width="30.63921959999999pt" height="24.65753399999998pt"/> the
curvature normal must point in the normal direction, the shape operator takes
as input a tangent vector and outputs another tangent vector (i.e., the change
in the unit normal must be tangential to the surface; no change can occur in
the normal direction itself).

Locally, the tangent vector space is two-dimensional spanned by basis vectors
<img alt="$\mathbf{e}_{1},\mathbf{e}_{2} \in  \mathbb{R}^{2}$" src="./tex/7caf22174becccaa731fcaa799d7340d.svg" align="middle" width="77.89937264999999pt" height="26.76175259999998pt"/> so we can think of the
shape operator as a mapping from <img alt="$\mathbb{R}^{2}$" src="./tex/3177e934cf575c08431076a1a5479ba5.svg" align="middle" width="18.424726649999986pt" height="26.76175259999998pt"/> to <img alt="$\mathbb{R}^{2}$" src="./tex/3177e934cf575c08431076a1a5479ba5.svg" align="middle" width="18.424726649999986pt" height="26.76175259999998pt"/>. As a differential operator,
the shape operator is a _linear operator_. This means we can represent its
action on a tangent vector <img alt="$\mathbf{v} = x_{1} \mathbf{e}_{1} + x_{2}\mathbf{e}_{2}$" src="./tex/7effb5a4c75b9fb50e8bcb8a655e1837.svg" align="middle" width="117.04311464999999pt" height="19.1781018pt"/>  as a matrix:

<p align="center"><img alt="$$&#10;S_\mathbf{p}(\mathbf{v}) = &#10;\left[&#10;\begin{array}{cc}&#10;S_\mathbf{p}(\mathbf{e}_{1})\cdot \mathbf{e}_{1} &amp; S_\mathbf{p}(\mathbf{e}_{1})\cdot \mathbf{e}_{2} \\&#10;S_\mathbf{p}(\mathbf{e}_{2})\cdot \mathbf{e}_{1} &amp; S_\mathbf{p}(\mathbf{e}_{2})\cdot \mathbf{e}_{2}&#10;\end{array}&#10;\right] \mathbf{v}&#10;$$" src="./tex/a847d0c122b10deb804b5cf91bf509e9.svg" align="middle" width="279.1205175pt" height="39.452455349999994pt"/></p>


Given <img alt="$\mathbf{r}_{1}$" src="./tex/61094f7c5090cb33cc734ef9b4d935ff.svg" align="middle" width="14.33791589999999pt" height="14.611878600000017pt"/> and <img alt="$\mathbf{r}_{2}$" src="./tex/aa13a17be39ba89a53599048c25da351.svg" align="middle" width="14.33791589999999pt" height="14.611878600000017pt"/> are the principal curvature directions (as unit 2D tangent
vectors) we can rotate our coordinate frame to align <img alt="$\mathbf{e}_{1}$" src="./tex/5e513183fa0be6505e7eb09ee857a2af.svg" align="middle" width="15.216900599999992pt" height="14.611878600000017pt"/> and <img alt="$\mathbf{e}_{2}$" src="./tex/d487296e678a76df356a18f1cf2e2bff.svg" align="middle" width="15.216900599999992pt" height="14.611878600000017pt"/> with the
principal curvature directions. The shape operator takes on a very special
form:

<p align="center"><img alt="$$&#10;S_\mathbf{p} = &#10;\left[\mathbf{r}_{1} \quad \mathbf{r}_{2}\right]&#10;\left[&#10;\begin{array}{cc}&#10;k_{1} &amp; 0 \\&#10;0 &amp; k^{2}&#10;\end{array}&#10;\right]&#10;\left[\mathbf{r}_{1} \quad \mathbf{r}_{2}\right]^{\mathsf T}&#10;$$" src="./tex/4fa4d113f75a8255f8e8a53ff5f582db.svg" align="middle" width="249.3057237pt" height="39.452455349999994pt"/></p>


> Consider why the off-diagonal terms are zero. Think about the _extremality_
> of the principal curvatures.

We have actually conducted an [eigen
decomposition](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix) on
the shape operator. Reading this progression backwards, the eigen decomposition
of the shape operator expressed in any basis will reveal:

 1. the principal curvatures as the eigen values, and
 2. the principal curvature directions as the eigen vectors.

<!--
Curvature for curves is the change in tangent vector under an arc-length. For a
given tangent direction on a surface, we extended this definition to define the
normal curvature as the curvature of the curve made by interesting surface with
a plane aligned with the chosen direction. Given two orthonormal tangent
directions <img alt="$\mathbf{u}$" src="./tex/129c5b884ff47d80be4d6261a476e9f1.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> and <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> (i.e., a local parameterization), let's collect the
normal curvature normal vectors:

<p align="center"><img alt="$$&#10;\begin{align}&#10;k_\mathbf{n}(\mathbf{u},\mathbf{p})\mathbf{n} &amp;= {\gamma}''_\mathbf{u}(\mathbf{p}) \\&#10;k_\mathbf{n}(\mathbf{v},\mathbf{p})\mathbf{n} &amp;= {\gamma}''_\mathbf{v}(\mathbf{p}).&#10;\end{align}&#10;$$" src="./tex/447488de72eefa6a34e2017819ade345.svg" align="middle" width="136.83982124999997pt" height="41.9471052pt"/></p>


If instead we equivalent consider the _change in normal vector_ for each sliced
curve, our curvature vectors will live in the orthogonal space: the tangent
space. 

<p align="center"><img alt="$$&#10;\begin{align}&#10;k_\mathbf{n}(\mathbf{u},\mathbf{p})\mathbf{v} &amp;= {\gamma}''_\mathbf{u}(\mathbf{p}) \\&#10;k_\mathbf{n}(\mathbf{v},\mathbf{p})\mathbf{u} &amp;= {\gamma}''_\mathbf{v}(\mathbf{p}).&#10;\end{align}&#10;$$" src="./tex/ec3872553101814de7f952fda4404847.svg" align="middle" width="136.5772782pt" height="41.9471052pt"/></p>


> Before we chose the normal direction by an angle <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/>, but for any tangent
> direction <img alt="$\mathbf{u}$" src="./tex/129c5b884ff47d80be4d6261a476e9f1.svg" align="middle" width="10.502226899999991pt" height="14.611878600000017pt"/> we can determine its correspond <img alt="${\varphi}$" src="./tex/c6badc6a64bc17f35f233989c2d6cbaa.svg" align="middle" width="10.75343279999999pt" height="14.15524440000002pt"/> so that <img alt="$k_\mathbf{n}({\varphi}) =&#10;&gt; k_\mathbf{n}(\mathbf{u})$" src="./tex/52a5082d1ea8ebcd574b2c73c828a677.svg" align="middle" width="116.81883674999999pt" height="24.65753399999998pt"/>.
-->

### Discrete curvatures on surfaces


#### Discrete mean curvature normal via discrete Laplace 

By now we are very familiar with the discrete Laplacian for triangle meshes:

<p align="center"><img alt="$$&#10;\Delta f \approx  \mathbf{M}^{-1} \mathbf{L} \mathbf{f},&#10;$$" src="./tex/2c98b11a27f5472ddf4cab80dc27dadd.svg" align="middle" width="104.53174214999999pt" height="17.399144399999997pt"/></p>

where <img alt="$\mathbf{M},\mathbf{L} \in \mathbb{R}^{n\times n}$" src="./tex/2fd2c03ca397d5c502202605ade1d029.svg" align="middle" width="95.11018605pt" height="26.17730939999998pt"/> are the mass and cotangent matrices respectively.

When applied to the vertex positions, this operator gives a point-wise (or
rather integral average) approximation of the mean curvature normal:

<p align="center"><img alt="$$&#10;H\mathbf{n} \approx  \mathbf{H} = \mathbf{M}^{-1} \mathbf{L} \mathbf{V} \in  \mathbb{R}^{n\times 3}.&#10;$$" src="./tex/c72dda0590d4586a4f51b74ca6e496a9.svg" align="middle" width="205.85359245pt" height="14.845497149999998pt"/></p>


Stripping the magnitude off the rows of the resulting matrix would give the
_unsigned_ mean curvature. To make sure that the sign is preserved we can check
whether each row in <img alt="$\mathbf{H}$" src="./tex/930b956ef51654e0669455a2cdd62fb5.svg" align="middle" width="14.794451099999991pt" height="22.55708729999998pt"/> agrees or disagrees with consistently oriented
per-vertex normals in <img alt="$\mathbf{N} \in  \mathbb{R}^{n\times 3}$" src="./tex/15de8d3756caf5b8d7e66cf1491fdf2e.svg" align="middle" width="71.71035959999999pt" height="26.76175259999998pt"/>. 

This connection between the Laplace operator and the mean curvature normal
provides additional understanding for its use as a geometric smoothing operator
(see "Computing Discrete Minimal Surfaces and Their Conjugates" [Pinkall and
Polthier 1993]).

#### Discrete Gaussian curvature via angle defect

On a discrete surface represented as a triangle mesh, curvature certainly can't
live on the flat faces. Moreover, Gaussian curvature can't live along edges
because we can always _develop_ the triangles on either side of an edge to the
plane without stretching them. In fact we can develop any arbitrarily long
chain of faces connected by edges so long as it doesn't form a loop or contain
all faces incident on a vertex. This hints that discrete Gaussian curvature
(like curvature for curves) must live at vertices.

Using the definition of Gaussian curvature in terms of the area on the Gauss
map in Equation <img alt="$(\ref{equ:gaussian-curvature-area})$" src="./tex/29719746ba75a77af37ed0f9b42180b1.svg" align="middle" width="30.63921959999999pt" height="24.65753399999998pt"/>: flat faces correspond
points on the Gauss map (contributing nothing), edges correspond to area-less
curves (traced by their [dihedral
angles](https://en.wikipedia.org/wiki/Dihedral_angle)), but vertices correspond
to spherical polygons connecting face normal-points. The area <img alt="${\Omega}$" src="./tex/9f531c9f3f1ebeef802ced46eabb0336.svg" align="middle" width="11.87217899999999pt" height="22.465723500000017pt"/> subtended on
the Gauss map is call the [solid
angle](https://en.wikipedia.org/wiki/Solid_angle). Conveniently, this area is
simply the [angle
defect](https://en.wikipedia.org/wiki/Angular_defect#Descartes.27_theorem) of
internal angles <img alt="${\theta}_f$" src="./tex/2d1a35f7ce73a328ee93ee9bfe9f0a8c.svg" align="middle" width="15.416760149999991pt" height="22.831056599999986pt"/> incident on the <img alt="$i$" src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/>-th vertex contributed by each <img alt="$f$" src="./tex/190083ef7a1625fbc75f243cffb9c96d.svg" align="middle" width="9.81741584999999pt" height="22.831056599999986pt"/>-th
incident face:

<p align="center"><img alt="$$&#10;{\Omega}_i = 2{\pi} - {\sum}\limits_{f \in  \text{faces(i)}} {\theta}_{if}.&#10;$$" src="./tex/723259007a16e34ff1b7e5cb445eb4f3.svg" align="middle" width="185.6342235pt" height="30.000197699999998pt"/></p>


!["Gaussian Curvature and Shell Structures" [Calladine
1986]](images/angle-defect.png)

Thus, our discrete analog of locally _integrated_ Gaussian curvature is given
as the angle defect at the <img alt="$i$" src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/>-th vertex. The local integral average (or
_pointwise_) discrete Gaussian curvature is the angle defect divided by the
local area associated with the <img alt="$i$" src="./tex/77a3b857d53fb44e33b53e4c8b68351a.svg" align="middle" width="5.663225699999989pt" height="21.68300969999999pt"/>-th vertex:

<p align="center"><img alt="$$&#10;K_i = \frac{2{\pi} - {\sum}\limits_{f \in  \text{faces(i)}} {\theta}_{if}}{A_i}.&#10;$$" src="./tex/17d205afb82b1e787ffb3a2a39d049ab.svg" align="middle" width="185.27574779999998pt" height="40.2896241pt"/></p>


By way of closing up the Gauss map, closed polyhedral surfaces (i.e., meshes)
will obey the
[Gauss-Bonnet](https://en.wikipedia.org/wiki/Gauss-Bonnet_theorem) in Equation 
<img alt="$(\ref{equ:gauss-bonnet})$" src="./tex/78ecd65ec74fbb72b055bdc8ecbc6383.svg" align="middle" width="30.63921959999999pt" height="24.65753399999998pt"/>, too:

<p align="center"><img alt="$$&#10;{\sum}\limits_{i=1}^n K_i = 2{\pi} {\chi}(\mathbf{S}).&#10;$$" src="./tex/0ead34ad413692d97b50e183f51ea550.svg" align="middle" width="143.53122299999998pt" height="28.013911200000003pt"/></p>


We can connect this to [Euler's
formula](https://en.wikipedia.org/wiki/Euler_characteristic) for polyhedra in our very first
assignment:

<p align="center"><img alt="$$&#10;\frac{1}{2{\pi}} {\sum}\limits_{i=1}^n K_i =  |V| - |E| + |F|,&#10;$$" src="./tex/8a021337d9eda79a54671df73fec95be.svg" align="middle" width="218.68865205pt" height="32.990165999999995pt"/></p>

where <img alt="$|V|, |E|, |F|$" src="./tex/1c4299bf1e3847240413930c7a0a782a.svg" align="middle" width="81.18724019999999pt" height="24.65753399999998pt"/> are the number of vertices, edges and faces respectively.



#### Approximation and eigen decomposition of the shape operator 

Alternatively, we can approximate all curvatures of a surface by locally
fitting an analytic surface and _reading_ off its curvature values. Since
planes have no curvature, the simplest type of analytic surface that will give
a non-trivial curvature value is a quadratic surface.

Thus, the algorithm proceeds as follows. For each vertex <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> of the given mesh,

 1.  gather a sampling of points in the vicinity. For simplicity, let's just
 grab all other vertices that share an edge with <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> or share an edge with a
 vertex that shares an edge with <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> (i.e., the "two-ring" of <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/>). For most
 sane meshes, this will provide enough points. Gather the positions of these
 <img alt="$k$" src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg" align="middle" width="9.075367949999992pt" height="22.831056599999986pt"/> points _relative_ to <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/> (i.e., <img alt="$\mathbf{v}_i - \mathbf{v}$" src="./tex/c24a23c5ffbaea48c8c7593a21285165.svg" align="middle" width="45.78079109999998pt" height="19.1781018pt"/>) into a matrix <img alt="$\mathbf{P} \in &#10; \mathbb{R}^{k\times 3}$" src="./tex/0c31361d39a47f2f1283726f5561e65f.svg" align="middle" width="68.97822194999999pt" height="27.91243950000002pt"/>.
 2. Next, we are going to define a quadratic surface as a height field above
 some two-dimensional plane passing through <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/>. Ideally, the plane is
 orthogonal to the normal at <img alt="$\mathbf{v}$" src="./tex/f6fc3ac36dff143d4aac9d145fadc77e.svg" align="middle" width="10.239687149999991pt" height="14.611878600000017pt"/>. To find such a plane, compute the
 [principal-component
 analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) of <img alt="$\mathbf{P}$" src="./tex/384591906555413c452c93e493b2d4ec.svg" align="middle" width="12.92230829999999pt" height="22.55708729999998pt"/>
 (i.e., conduct eigen decomposition on <img alt="$\mathbf{P}^{\mathsf T} \mathbf{P}$" src="./tex/e271558c00da98acfd7fa2a6d27d49f5.svg" align="middle" width="35.01712169999999pt" height="27.91243950000002pt"/>). Let <img alt="$\mathbf{S} \in  \mathbb{R}^{k \times &#10; 2}$" src="./tex/4a33bc01ee03345e0a7a058b12680ff0.svg" align="middle" width="66.55814219999998pt" height="27.91243950000002pt"/> be the coefficients for two most principal directions (call them the <img alt="$u$" src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg" align="middle" width="9.41027339999999pt" height="14.15524440000002pt"/>-
 and <img alt="$v$" src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg" align="middle" width="8.55786029999999pt" height="14.15524440000002pt"/>- directions) corresponding to each point in <img alt="$\mathbf{P}$" src="./tex/384591906555413c452c93e493b2d4ec.svg" align="middle" width="12.92230829999999pt" height="22.55708729999998pt"/>, and let <img alt="$\mathbf{B} \in &#10; \mathbb{R}^{k}$" src="./tex/ce6e2d702a9c522b607b1a1add2fb6cd.svg" align="middle" width="52.67676644999999pt" height="27.91243950000002pt"/> be the "height" of each point in the least principal direction (call
 it the <img alt="$w$" src="./tex/31fae8b8b78ebe01cbfbe2fe53832624.svg" align="middle" width="12.210846449999991pt" height="14.15524440000002pt"/>-direction).
 3. An quadratic function as a height-field surface passing through the origin
 is given by:
<p align="center"><img alt="$$&#10;w = a_{1} u + a_{2} v + a_{3} u^{2} + a_{4}uv + a_{5}v^{2}.&#10;$$" src="./tex/db56e821d2b5b5510ecaa175bb4ac666.svg" align="middle" width="268.0307982pt" height="16.66852275pt"/></p>

 We have <img alt="$k$" src="./tex/63bb9849783d01d91403bc9a5fea12a2.svg" align="middle" width="9.075367949999992pt" height="22.831056599999986pt"/> sets of <img alt="$u,v$" src="./tex/cfecde842a36413fb233cf4913fbcb8f.svg" align="middle" width="25.27401689999999pt" height="14.15524440000002pt"/> values and <img alt="$w$" src="./tex/31fae8b8b78ebe01cbfbe2fe53832624.svg" align="middle" width="12.210846449999991pt" height="14.15524440000002pt"/> values. Treat this as a
 least-squares fitting problem and solve for the 5 unknown coefficients.
 (`igl::pinv` is good for solving this robustly).
 4. Each element of the shape operator for the graph of a quadratic function
 over the plane has a closed form expression. You need to derive these by hand.
 Just kidding. The shape operator can be constructed as the product of two
 matrices:
<p align="center"><img alt="$$&#10;S = -&#10;\left[&#10;\begin{array}{cc}&#10;e &amp; f \\&#10;f &amp; g&#10;\end{array}&#10;\right]&#10;\left[&#10;\begin{array}{cc}&#10;E &amp; F \\&#10;F &amp; G&#10;\end{array}&#10;\right]^{-1}&#10;$$" src="./tex/c300d192ba79563a6b074472d88063c8.svg" align="middle" width="214.13473125pt" height="42.80407395pt"/></p>

 known as the second and first fundamental forms respectively. The entries of
 these matrices categorize the stretch and bending in each direction:
<p align="center"><img alt="$$&#10;E = 1+a_{1}^{2} \\&#10;F = a_{1}a_{2}  \\&#10;G = 1+a_{2}^{2} \\&#10;e = \frac{2a_{3}}{\sqrt{a_{1}^{2} + 1 + a_{2}^{2}}} \\&#10;f = \frac{a_{4}}{\sqrt{a_{1}^{2} + 1 + a_{2}^{2}}} \\&#10;g = \frac{2a_{5}}{\sqrt{a_{1}^{2} + 1 + a_{2}^{2}}}&#10;$$" src="./tex/037966b300e7c24702df8e3528cd817a.svg" align="middle" width="617.9080941pt" height="40.289634pt"/></p>

 See Table 1 of "Estimating Differential Quantities Using Polynomial Fitting of
 Osculating Jets" [Cazals & Pouget 2003] to double check for typos :-).
 5. Eigen decomposition of <img alt="$S$" src="./tex/e257acd1ccbe7fcb654708f1a866bfe9.svg" align="middle" width="11.027402099999989pt" height="22.465723500000017pt"/> reveals the principal curvatures <img alt="$k_{1}$" src="./tex/df2a706e439b4632680202a540f27fd0.svg" align="middle" width="15.11042279999999pt" height="22.831056599999986pt"/> and <img alt="$k_{2}$" src="./tex/bce145d2ce5a01a17e36278910bfa8a4.svg" align="middle" width="15.11042279999999pt" height="22.831056599999986pt"/>
 _and_ the principal tangent directions (in the <img alt="$uv$" src="./tex/ee5f11272c9cd93256bbf7ba019c3953.svg" align="middle" width="17.96813369999999pt" height="14.15524440000002pt"/> PCA basis).
 6. Lift the principal tangent directions back to world <img alt="$\mathbb{R}^{3}$" src="./tex/fabbfebc4049d77e28eefb36851e7538.svg" align="middle" width="18.424726649999986pt" height="26.76175259999998pt"/> coordinates.

## Tasks

[Download](https://archive.org/details/ElementaryDifferentialGeometry) Barret
O'Neill's book. This is my go-to differential geometry book. The section on
curvature and the shape operator should help resolve questions and fill in
missing proofs above.

### Blacklist

 - `igl::gaussian_curvature`
 - `igl::internal_angles` (or any of the other overloads)
 - `igl::principal_curvatures`

### Whitelist

 - `igl::adjacency_matrix.h`
 - `igl::cotmatrix`
 - `igl::invert_diag`
 - `igl::massmatrix`
 - `igl::per_vertex_normals`
 - `igl::pinv`
 - `igl::slice`
 - `igl::sort`
 - `igl::squared_edge_lengths`

### `src/mean_curvature.cpp`
Compute the discrete mean curvature at each vertex of a mesh (`V`,`F`) by
taking the signed magnitude of the mean curvature normal as a _pointwise_ (or
_integral average_) quantity.

### `src/internal_angles.cpp`
Given (squared) edge-lengths of a triangle mesh `l_sqr` compute the internal
angles at each corner (a.k.a. wedge) of the mesh.

### `src/angle_defect.cpp`
Compute the discrete angle defect at each vertex of a triangle mesh
(`V`,`F`), that is, the _locally integrated_ discrete Gaussian
curvature.

### `src/principal_curvatures.cpp`
Approximate principal curvature values and directions locally by considering
the two-ring neighborhood of each vertex in the mesh (`V`,`F`).
