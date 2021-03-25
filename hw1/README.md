# Discrete and Smooth Curves


### Problem 1
In this problem, we introduce continuous and discrete methods in variational calculus, one of the main tools of the differential geometry toolbox.

**a)** Suppose you are given a regular plane curve $\gamma: [0,1] \to \mathbb{R}^2$, and take $\mathbf{v}: [0,1] \to \mathbb{R}^2$ to be a vector field along $\gamma$. Recall that the arc length of $\gamma$ is given by
```math
s[\gamma] = \int_0^1 \lVert\gamma'(t)\rVert_2\ dt.
```
We can think of $\gamma_h(t) := \gamma(t) + h\mathbf{v}(t)$ to be a displacement of $\gamma$ along $\mathbf{v}$. Differentiate $f(h) := s[\gamma + h\mathbf{v}]$ with respect to $h$ at $h=0$ to yield an expression for $\left. \frac{d}{dh} s[\gamma + h\mathbf{v}] \right|_{h=0}$.

**answer:**
```math
\begin{align*}
  \left. \frac{d}{dh} s[\gamma + h\mathbf{v}] \right|_{h=0}
    &= \left. \frac{d}{dh} \left( \int_0^1 \lVert\gamma'(t) + h\mathbf{v}'(t)\rVert\ dt \right) \right|_{h=0}\\
    &= \int_0^1 \left. \tfrac{d}{dh} \lVert\gamma'(t) + h\mathbf{v}'(t)\rVert \right|_{h=0}\ dt \\
    &= \int_0^1 \left. \frac{\gamma'(t) \cdot \mathbf{v}'(t) + h \lVert\mathbf{v}'(t)\rVert}{\lVert\gamma'(t) + h\mathbf{v}'(t)\rVert} \right|_{h=0}\ dt \\
    &= \int_0^1 \frac{\gamma'(t)}{\lVert\gamma'(t)\rVert} \cdot \mathbf{v}'(t)\ dt
\end{align*}
```

**b)** You can think of each $\mathbf{v}$ as an infinitesimal displacement (a "variation") of the entire curve $\gamma$ at once. Explain how the derivative you took in part **a** can be thought of as a directional derivative of arc length in the $\mathbf{v}$ "direction." In variational calculus, this derivative is known as the Gateaux or variational derivative of $s[\cdot]$.

**answer:** First, note that we have:
```math
\begin{align*}
  \left. \frac{d}{dh} s[\gamma + h\mathbf{v}] \right|_{h=0}
    &= \left. \frac{d}{dh} f(h) \right|_{h=0} \\
    &= f'(0) \\
    &= \lim_{h \to 0} \frac{f(h) - f(0)}{h} \\
    &= \lim_{h \to 0} \frac{s[\gamma + h\mathbf{v}] - s[\gamma]}{h}.
\end{align*}
```
We can think of the arc length $s$ as the function of smooth (or at least twice differentiable) functions defined on $[0,1]$ i.e. $s: \mathcal{C}^{\infty}\left([0,1],\mathbb{R}^2\right) \to \mathbb{R}$. Then, we can compute its variational derivative at $\gamma$ in the direction $\mathbf{v}$:
```math
ds[\gamma; \mathbf{v}] = \lim_{h \to 0} \frac{s[\gamma + h\mathbf{v}] - s[\gamma]}{h}.
```

**c)** Suppose $\mathbf{v}(0) = \mathbf{v}(1) = 0$. Define a vector-valued function $\mathbf{w}(s)$ so that
```math
\left. \frac{d}{dh} s[\gamma + h\mathbf{v}] \right|_{h=0} = \int_0^{s(1)} \mathbf{v}(s^{-1}(\bar{s})) \cdot \mathbf{w}(\bar{s})\ d\bar{s},
```
where $s(t) = \int_0^t \lVert\gamma'(\bar{t})\rVert_2\ d\bar{t}$ on the right-hand side is the arc length function and $\mathbf{w}$ can be written in terms of the curvature and Frenet frame of $\gamma$.

**answer:** By $\mathbf{t}$ and $\mathbf{n}$ let us denote the tangent and normal vectors respectively of the Frenet frame of $\gamma$. Also let $\kappa(t)$ be a curvature of $\gamma$ at the point $t$. Now, let us define a vector-valued function $\mathbf{w}$ as $\mathbf{w}(t) = -\kappa(t) \mathbf{n}(t)$. Then, we have
```math
\begin{align*}
  \left. \frac{d}{dh} s[\gamma + h\mathbf{v}] \right|_{h=0}
    &= \int_0^1 \mathbf{v}'(t) \cdot \tfrac{\gamma'(t)}{\lVert\gamma'(t)\rVert}\ dt \\
    &= \left. \mathbf{v}(t) \cdot \tfrac{\gamma'(t)}{\lVert\gamma'(t)\rVert} \right|_0^1 - \int_0^1 \mathbf{v}(t) \cdot \left(\tfrac{d}{dt} \tfrac{\gamma'(t)}{\lVert\gamma'(t)\rVert}\right)\ dt \\
    &= \int_0^1 \mathbf{v}(t) \cdot \left(-\tfrac{d}{dt} \tfrac{\gamma'(t)}{\lVert\gamma'(t)\rVert}\right)\ dt \\
    &= \int_0^{s(1)} \mathbf{v}(s^{-1}(\bar{s})) \cdot \left(-\tfrac{d}{d\bar{s}} \mathbf{t}(\bar{s})\right)\ d\bar{s} \\
    &= \int_0^{s(1)} \mathbf{v}(s^{-1}(\bar{s})) \cdot \left(-\kappa(\bar{s}) \mathbf{n}(\bar{s})\right)\ d\bar{s} \\
    &= \int_0^{s(1)} \mathbf{v}(s^{-1}(\bar{s})) \cdot \mathbf{w}(\bar{s})\ d\bar{s}.
\end{align*}
```


### Problem 2
In this problem, you will develop a notion of discrete curvature of a plane curve.

**a)** Suppose we have a discrete curve given by a series of points $\mathbf{x}_1,\ldots,\mathbf{x}_n \in \mathbb{R}^2$. You can think of the vertex positions as parameterized by a vector $\mathbf{x} \in \mathbb{R}^{2n}$. Define an arc length functional $s(\mathbf{x}): \mathbb{R}^{2n} \to \mathbb{R}_{+}$; for convenience, it is acceptable to notate $s(\mathbf{x}) = s(\mathbf{x}_1,\ldots,\mathbf{x}_n)$.

**answer:** Let us define the arc length functional $s$ as a sum of all linear segments between the points $\mathbf{x}_1,\ldots,\mathbf{x}_n$:
```math
s(\mathbf{x}) = \sum_{i=1}^{n-1} \lVert \mathbf{x}_{i+1} - \mathbf{x}_{i} \rVert
```

**b)** Suppose $1 < i < n$. Write an expression for the gradient $\nabla_{\mathbf{x}_i} s$ of $s(\cdot)$ with respect to $\mathbf{x}_i$ and show that its norm is $2\sin\frac{\theta}{2}$, where $\theta$ is the turning angle between the two segments adjacent to $\mathbf{x}_i$.

**answer:**
```math
\nabla_{\mathbf{x}_i} s(\mathbf{x}) = \frac{\mathbf{x}_{i} - \mathbf{x}_{i-1}}{\lVert \mathbf{x}_{i} - \mathbf{x}_{i-1} \rVert} - \frac{\mathbf{x}_{i+1} - \mathbf{x}_{i}}{\lVert \mathbf{x}_{i+1} - \mathbf{x}_{i} \rVert}
```
Now, let us denote by $\mathbf{t}_{i}$ the unit tangent vector to the linear segment between the points $\mathbf{x}_{i}$ and $\mathbf{x}_{i+1}$ for $i=1,2,\ldots,n-1$, i.e:
```math
\mathbf{t}_{i} = \frac{\mathbf{x}_{i+1} - \mathbf{x}_{i}}{\lVert \mathbf{x}_{i+1} - \mathbf{x}_{i} \rVert}.
```
Then, we have $\nabla_{\mathbf{x}_i} s(\mathbf{x}) = \mathbf{t}_{i-1} - \mathbf{t}_{i}$ and
```math
\begin{align*}
  \rVert \nabla_{\mathbf{x}_i} s(\mathbf{x}) \lVert^2
    &= \rVert \mathbf{t}_{i-1} - \mathbf{t}_{i} \lVert^2 \\
    &= \rVert \mathbf{t}_{i-1} \lVert^2 + \rVert \mathbf{t}_{i}\lVert^2 - 2\, \mathbf{t}_{i-1} \cdot \mathbf{t}_{i} \\
    &= 2\left(1 - \mathbf{t}_{i-1} \cdot \mathbf{t}_{i} \right) \\
    &= 2\left(1 - \cos\theta \right) \\
    &= 2\left(\left(\cos^2\tfrac{\theta}{2} + \sin^2\tfrac{\theta}{2}\right) - \left(\cos^2\tfrac{\theta}{2} - \sin^2\tfrac{\theta}{2}\right) \right) \\
    &= 4\sin^2\tfrac{\theta}{2},
\end{align*}
```
which proves that the norm of $\nabla_{\mathbf{x}_i} s$ equals $2\sin\tfrac{\theta}{2}$.

**c)** Take a look at `discreteCurve.jl`.
The curve generates an $n \times 2$ array representing $n$ points on a discrete two-dimensional curve. Modify the code to plot the derivative you computed in part **b**.

**answer:** See `discreteCurve.jl`.

**d)** Propose a measure of discrete (unsigned) per-vertex curvature of a 2D discrete curve based on your answers to **1.c** and **2.b**, and draw the curve colored by this value. For this, fill in code to compute $\kappa$.

**answer:** See `discreteCurve.jl`.

**e)** For sufficiently small $h > 0$, one simple way to decrease the length of the curve would be to replace each point $\mathbf{x}_i$ with a new point $\mathbf{x}_i' := \mathbf{x}_i − (\nabla_{\mathbf{x}_i} s) h$, where $\nabla_{\mathbf{x}_i} s$ is the derivative of $s$ with respect to $\mathbf{x}_i$.
Implement this forward integration scheme with the endpoints fixed, and make sure that if you iterate enough times the curve approximates a straight line. What happens if $h$ is too large?

**answer:** See `discreteCurve.jl`.
For small $h > 0$ applying the gradient descend scheme to the curve with fixed endpoints $\mathbf{x}_1$ and $\mathbf{x}_n$ leads to the fully straigtened curve when run enough times. Too small values of $h$ results in smaller convergence rate and longer optimization times. On the other hand, too large values of $h$ can prevent the algorithm to converge.


### Problem 3
In this problem, you will implement part of the "Discrete Elastic Rods" paper.

**a)** Add code to compute the $(n−2) \times 3$ array binormal, which contains the Darboux vector $(\kappa\mathbf{b})_i$ for each vertex $i$ except the first and last.

**answer:** See `bishopFrame.jl`.

**b)** We provide simple code for animating different initial choices of the Bishop frame $(\mathbf{u},\mathbf{v},\mathbf{t})$ on the first segment. Add code to fill in $\mathbf{u}$ and $\mathbf{v}$ along the rest of the curve.

**answer:** See `bishopFrame.jl`.

**c)** In the final part of the code, we prescribe a material frame on the curve. Add code to compute the twist energy and its gradient. Perform gradient descent on the angles between the material frame and the Bishop frame and observe how the material frame "untwists." Make sure that if you run your code long enough the material frame aligns with the natural Bishop frame.

**answer:** See `bishopFrame.jl`.


### References
- Differential Geometry of Curves and Surfaces, do Carmo, Courier Dover Publications, 2016
- Curves and Surfaces, Montiel and Ros, AMS, 2009
- Elasticity and Geometry, Audoly and Pomeau, Oxford University Press, 2010
- There is more than one way to frame a curve, Bishop, The American Mathematical Monthly 82(3), 1975
- The Resultant Electric Moment of Complex Molecules, Eyring, Physical Review 39(4), 1932
- The Curvatures of Regular Curves and Euclidean Invariants of their Derivatives, Gutkin, arXiv 1007.2960, 2010
- Discrete Elastic Rods, Bergou et al., SIGGRAPH, 2008
- Discrete Viscous Threads, Bergou et al., SIGGRAPH, 2010
