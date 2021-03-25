# Surfaces and Curvature


### Problem 1
In this problem, you'll do a bit of calculus to see how the operators we talk about in class work on a simple manifold.
Consider the map $f: \mathbb{R}^2 \to \mathbb{R}^3$ defined by
```math
f(u,v) = \frac{1}{u^2+v^2+1}(2u,2v,u^2+v^2-1).
```

**a)** Verify that for all $(u,v) \in \mathbb{R}^2$, $f(u,v)$ lies on the unit sphere.

**answer:** To show that $f(u,v)$ lies on the unit sphere we prove that $\lVert f(u,v) \rVert = 1$:
```math
\begin{align*}
  \lVert f(u,v) \rVert
    &= \frac{\lVert(2u,2v,u^2+v^2-1)\rVert}{u^2+v^2+1} \\
    &= \frac{\sqrt{u^4+v^4+1+2u^2v^2+2u^2+2v^2}}{u^2+v^2+1} \\
    &= \frac{\sqrt{(u^2+v^2+1)^2}}{u^2+v^2+1} \\
    &= 1.
\end{align*}
```

**b)** Let $\mathbf{p} = (u_0,v_0)$ be a point in $\mathbb{R}^2$, and let $\gamma: (-\varepsilon,\varepsilon) \to \mathbb{R}^2$ be a curve with $\gamma(0) = \mathbf{p}$, and $\gamma'(0) = \mathbf{v}$.
Recall that the differential of a function $f$ at a point $\mathbf{p}$ is a linear map $df_{\mathbf{p}}(\mathbf{v}) = (f\circ \gamma)'(0)$.
Compute the differential map $df_{\mathbf{p}}$ at $\mathbf{w} = \mathbf{e}_1 w^1 + \mathbf{e}_2 w^2$.

**answer:** From the definition of a differential:
```math
\begin{align*}
  df_{\mathbf{p}}(\mathbf{w})
    &= df_{\mathbf{p}}(\mathbf{e}_1 w^1 + \mathbf{e}_2 w^2) \\
    &= w^1 df_{\mathbf{p}}(\mathbf{e}_1) +
       w^2 df_{\mathbf{p}}(\mathbf{e}_2) \\
    &= w^1 \tfrac{\partial f}{\partial u}(\mathbf{p}) +
       w^2 \tfrac{\partial f}{\partial v}(\mathbf{p}).
\end{align*}
```
Thus:
```math
df_{\mathbf{p}}(\mathbf{w})
  = \frac{2w^1}{(u_0^2+v_0^2+1)^2} \begin{pmatrix} -u_0^2+v_0^2+1 \\ -2u_0v_0 \\ 2u_0 \end{pmatrix} +
    \frac{2w^2}{(u_0^2+v_0^2+1)^2} \begin{pmatrix} -2u_0v_0 \\ u_0^2-v_0^2+1 \\ 2v_0 \end{pmatrix}.
```

**c)** Recall that the Gauss map of a surface $\mathcal{M}$ is a function $n: \mathcal{M} \to S^2$. Given a parameterization of a surface, there is a simple way to obtain the Gauss map $n(\mathbf{p})$. What is the Gauss map induced by $f$?

**answer:** Gauss map of a oriented surface given the parametrization $f$ can be defined using cross product:
```math
n(\mathbf{p}) =
\frac{\tfrac{\partial f}{\partial u}(\mathbf{p}) \times \tfrac{\partial f}{\partial v}(\mathbf{p})}
{\lVert \tfrac{\partial f}{\partial u}(\mathbf{p}) \times \tfrac{\partial f}{\partial v}(\mathbf{p}) \rVert},
```
which is equal to $-f(\mathbf{p})$.

**d)** Compute the differential $dn_{\mathbf{p}}(\mathbf{v})$ of the Gauss map.

**answer:** From the previous part we know that $n(\mathbf{p}) = -f(\mathbf{p})$. Therefore:
```math
dn_{\mathbf{p}}(\mathbf{v})
  = \frac{-2v^1}{(u_0^2+v_0^2+1)^2} \begin{pmatrix} -u_0^2+v_0^2+1 \\ -2u_0v_0 \\ 2u_0 \end{pmatrix} +
    \frac{-2v^2}{(u_0^2+v_0^2+1)^2} \begin{pmatrix} -2u_0v_0 \\ u_0^2-v_0^2+1 \\ 2v_0 \end{pmatrix}.
```


### Problem 2
Recall the Taubin matrix defined in class for approximating mesh curvature
```math
M_{\mathbf{p}} = \frac{1}{2\pi} \int_{-\pi}^{\pi} \kappa_\theta \mathbf{t}_\theta \mathbf{t}_\theta^\top\ d\theta.
```

**a)** Prove that the surface normal at $\mathbf{p}$ is an eigenvector of $M_{\mathbf{p}}$. What is the corresponding eigenvalue?

**answer:** Let us compute $M_{\mathbf{p}}\mathbf{n}$:
```math
\begin{align*}
  M_{\mathbf{p}}\mathbf{n}
    &= \frac{1}{2\pi} \int_{-\pi}^{\pi} \kappa_\theta \mathbf{t}_\theta \mathbf{t}_\theta^\top \mathbf{n}\ d\theta \\
    &= \frac{1}{2\pi} \int_{-\pi}^{\pi} \kappa_\theta \mathbf{t}_\theta (\mathbf{t}_\theta \cdot \mathbf{n})\ d\theta \\
    &= \frac{1}{2\pi} \int_{-\pi}^{\pi} 0\ d\theta \\
    &= 0.
\end{align*}
```
This shows that the surface normal vector $\mathbf{n}$ is in the kernel of $M_{\mathbf{p}}$ i.e. it is its eigenvector with corresponding eigenvalue $0$.

**b)** Show that the other two eigenvectors are the principal curvature directions. What are the corresponding eigenvalues?

**answer:** First, note that $\mathbf{t}_1^\top M_{\mathbf{p}} \mathbf{t}_2 = \mathbf{t}_2^\top M_{\mathbf{p}} \mathbf{t}_1$ and let us compute this value:
```math
\begin{align*}
  \mathbf{t}_1^\top M_{\mathbf{p}} \mathbf{t}_2
    &= \frac{1}{2\pi} \int_{-\pi}^{\pi} \kappa_\theta \cos\theta\sin\theta\ d\theta \\
    &= \frac{\kappa_1}{2\pi} \int_{-\pi}^{\pi} \cos^3\theta\sin\theta\ d\theta +
       \frac{\kappa_2}{2\pi} \int_{-\pi}^{\pi} \cos\theta\sin^3\theta\ d\theta \\
    &= 0.
\end{align*}
```
Recall that the principal curvature directions $\mathbf{t}_1, \mathbf{t}_2$ make an orthonolmal basis of the tangent plane. From the part **a** we know the normal vector $\mathbf{n}$ belongs to the kernel of $M_{\mathbf{p}}$. These facts together with the result above prove that the vectors $\mathbf{t}_1, \mathbf{t}_2$ are eigenvectors of $M_{\mathbf{p}}$. Now, we can find its corresponding eigenvalues $\lambda_1, \lambda_2$ with the following calculations:
```math
\begin{align*}
  \lambda_1
    &= \mathbf{t}_1^\top M_{\mathbf{p}} \mathbf{t}_1 \\
    &= \frac{1}{2\pi} \int_{-\pi}^{\pi} \kappa_\theta \cos^2\theta\ d\theta \\
    &= \frac{\kappa_1}{2\pi} \int_{-\pi}^{\pi} \cos^4\theta\ d\theta +
       \frac{\kappa_2}{2\pi} \int_{-\pi}^{\pi} \cos^2\theta\sin^2\theta\ d\theta \\
    &= 3 \frac{\kappa_1}{8} + \frac{\kappa_2}{8} \\
  \lambda_2
    &= \mathbf{t}_2^\top M_{\mathbf{p}} \mathbf{t}_2 \\
    &= \frac{1}{2\pi} \int_{-\pi}^{\pi} \kappa_\theta \sin^2\theta\ d\theta \\
    &= \frac{\kappa_1}{2\pi} \int_{-\pi}^{\pi} \cos^2\theta\sin^2\theta\ d\theta +
       \frac{\kappa_2}{2\pi} \int_{-\pi}^{\pi} \sin^4\theta\ d\theta \\
    &= \frac{\kappa_1}{8} + 3 \frac{\kappa_2}{8} \\
\end{align*}
```


### Problem 3
In this problem, you will compute and display (pointwise, not integrated) Gaussian curvature on a triangle mesh. Since there are many approximations for discrete Gaussian curvature, choose any two and fill in the relevant portions of `gaussCurvature.jl`.
In one or two sentences, compare your two choices of curvature. Are there situations in which they behave differently?

**answer:** See `gaussCurvature.jl`.
Implemented three approaches:
- `gaussianCurvature1`: Gaussian curvature estimation by structure preservation (Gauss-Bonnet theorem),
- `gaussianCurvature2`: Gaussian curvature estimation by the Taubin matrix,
- `gaussianCurvature3`: Gaussian curvature estimation by approximation of the second fundamental form.


### Problem 4
In this problem, you will develop a notion of pointwise mean curvature ona triangle mesh. Take a look at `meanCurve.jl` for starter code.

**a)** Complete the function `surfaceArea` which computes the surface area of a triangle mesh from the vertices and triangles.

**answer:** See `meanCurve.jl`.

**b)** Complete the function `cotLaplacian` that computes a sparse matrix $L$ such that $\nabla_{\mathbf{p}}A = \frac{1}{2} L \cdot \mathbf{p} \in \mathbb{R}^{|V| \times 3}$, where $A$ is the surface area from the previous part, $\mathbf{p} \in \mathbb{R}^{|V| \times 3}$ contains vertex positions, and $L \in \mathbb{R}^{|V|\times|V|}$ depends on $\mathbf{p}$ and the topology of the mesh.

**answer:** See `meanCurve.jl`.

**c)** Suppose we want to check that our cotangent Laplacian is indeed approximating the gradient of area. One way is to compute $\nabla_{\mathbf{p}}A$ via divided differences on $A$ with respect to point positions. Complete the function `dividedDifferences` and show that the error between this approximation and the true gradient you computed in part **b** is small.

**answer:** See `meanCurve.jl`.

**d)** The barycentric area associated to a vertex is $1/3$ times the sum of triangle areas adjacent to that vertex. Complete the function `barycentricArea`, and argue that the sum of barycentric areas over all vertices is the surface area.

**answer:** See `meanCurve.jl`.

**e)** Combine code from the previous parts to approximate a per-vertex pointwise mean curvature on the mesh. Fill in the `meanCurvature` function.

**answer:** See `meanCurve.jl`.


### Problem 5
Now, you will use code from problem **4** to implement the mean curvature flow algorithm described in "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow." See the relevant portions of `meanCurve.jl` for starter code.

**a)** Take $M \in \mathbb{R}^{|V| \times |V|}$ to be a diagonal matrix of barycentric areas from problem **4.c**. Notice that $M$ and $L$ are functions of $\mathbf{p}$: $M, L : \mathbb{R}^{|V| \times 3} \to \mathbb{R}^{|V| \times |V|}$. Based on our discussion of the mean curvature normal, how do you expect the following ODE to evolve $\mathbf{p}$ in time $t \geq 0$:
```math
\frac{d\mathbf{p}}{dt} = -M(\mathbf{p})^{−1} \cdot L(\mathbf{p}) \cdot \mathbf{p}.
```

**answer:** We know that $\triangle\mathbf{p} = 2H\mathbf{n}$. On the other hand, the latter expression is equal to the variational derivative of the surface area. Thus, the ODE needs to smooth $\mathbf{p}$ in time.

**b)** Suppose we wish to approximate $\mathbf{p}(t+\tau)$ given $\mathbf{p}(t)$. One simple way is to solve the following divided difference approximation for $\mathbf{p}(t+\tau)$:
```math
\frac{\mathbf{p}(t+\tau)-\mathbf{p}(t)}{\tau} \approx −M(\mathbf{p}(t))^{−1} \cdot L(\mathbf{p}(t)) \cdot \mathbf{p}(t).
```
Implement this approximation in `curvatureFlowEuler`. What happens if $\tau$ is too large?

**answer:** See `meanCurve.jl`.
Too large values of $\tau$ give poor approximations and often prevent the algorithm from convergence.

**c)** An alternative (semi-)implicit integrator uses a different approximation:
```math
\frac{\mathbf{p}(t+\tau)−\mathbf{p}(t)}{\tau} \approx −M(\mathbf{p}(t))^{−1} \cdot L(\mathbf{p}(t)) \cdot \mathbf{p}(t+\tau).
```
Implement this approximation in `curvatureFlowImplicit`. What happens if $\tau$ is too large?

**answer:** See `meanCurve.jl`.
Although a semi-implicit integrator is more stable it is still sensitive to large values of $\tau$. Using too large $\tau$ values can prevent the algorithm from convergence.

**d)** A fully-implicit integrator would use the following approximation:
```math
\frac{\mathbf{p}(t+\tau)−\mathbf{p}(t)}{\tau} \approx −M(\mathbf{p}(t+\tau))^{-1} \cdot L(\mathbf{p}(t+\tau)) \cdot \mathbf{p}(t+\tau).
```
Speculate in words why this formula is not used as often as the previous two.

**answer:** The explicit and semi-implicit formulas yield linear equations on $\mathbf{p}(t+\tau)$. This is, however, no longer true for the fully-implicit integrator. Due to the computational cost of solving non-linear equations it is not used as often as the previous two methods.


### Problem 6
In the previous problem, you might have seen that mean curvature flow produces singularities and sharp features when it is run for too long. The paper "Can Mean-Curvature Flow be Modified to be Non-Singular?" (Kazhdan, Solomon, and Ben-Chen) proposes a solution to this problem. Implement their approach and report on the results.

**answer:** See `meanCurve.jl`.
The solution to the problem of singularities is given in the section "3.2 Conformalizing the Metric" of the paper, mainly:
> For the finite-elements discretization, the implementation of this modification is trivial. Instead of requiring that the stiffness matrix be computed anew at each time-step (as is required by traditional mean-curvature flow), the modified flow simply reuses the stiffness matrix from time $t = 0$.


### References
- Differential Geometry of Curves and Surfaces, do Carmo, Courier Dover Publications, 2016
- Curves and Surfaces, Montiel and Ros, AMS, 2009
- Curvatures of Smooth and Discrete Surfaces, Sullivan, arXiv 0710.4497, 2007
- Discrete Differential-Geometry Operators for Triangulated 2-Manifolds, Meyer et al., VisMath, 2002
- Estimating the Tensor of Curvature of a Surface from a Polyhedral Approximation, Taubin, ICCV, 1995
- Estimating curvatures and their derivatives on triangle meshes, Rusinkiewicz, 3DPVT, 2004
- Implicit fairing of irregular meshes using diffusion and curvature flow, Desbrun et al., SIGGRAPH, 1999
- Can Mean-Curvature Flow be Modified to be Non-Singular?, Kazhdan et al., arXiv 1203.6819, 2012
