
---
title: 'Optimization on Stiefel manifolds'
date: 2020-08-24
permalink: /posts/2020/08/blog-post-1/
tags:
  - optimization
  - differential topology
  - differential geometry
---

Here I present an interesting way of formulating the constraints in an optimization problem, namely in terms of differential geometry.

Note: This post is based on the great notes I found [http://noodle.med.yale.edu/~hdtag/notes/steifel_notes.pdf](here). My goal is to provide an intuition of the problem, a high-level description of the solution as well as a neat strategy of how to formulate and approach constrained optimization problems. If you want to go deeper into the math you should check the notes. 

Problem formulation
-----

Suppose we have following optimization problem:
$$\displaystyle \min _{X \in \mathcal{R}^{n \times p}} F(X) \text { s.t. } X^{T} X=I$$ where $F : \mathcal{R}^{n \times p} \to \mathcal{R}$ is a score function.

Here, we may be interested in the matrix $X$, or we may want to find $p$ orthogonal vectors which are optimal wrt. to a score function. This problem has applications for example in RNNs.

We are going to use a gradient-based approach on the manifold spanned by the constraint. Why not use gradient descent directly in $\mathcal{R}^{n \times p}$? Well, there are many reasons, for example curse of dimensionality of the convergence problems in high dimensions. We can circumvent these issues by working with a manifold instead.  

The first step is to show that the space of orthogonal matrices is a manifold. This is a very theoretical and tedious step and, unless you are a pure mathematician with experience in the field, you should check if someone has already found that for you. If not, you should try another approach. In our case, a proof is given in <cite>[1]</cite>, Example 6.10. This manifold is called a <i>Stiefel manifold</i> and is a submanifold of $\mathbb{R}^{n \times p}$.

The immediate thing to do is to identify the tangent space, i.e. the space of vectors which are orthogonal to an element in our manifold. We find that each tangent vector satisfies $Z^{T} X+X^{T} Z=0$.

Remember: Differential topology is very abstract and mathematicians often work with abstract objects. But in order to do computations on a manifold, we need to find a representation of these objects. On the other hand, this comes in handy when we want functions which are invariant under reparametrization, such as ([https://en.wikipedia.org/wiki/Jeffreys_prior](https://en.wikipedia.org/wiki/Jeffreys_prior))[Jeffreys prior], which is uniform on the statistical manifold. 

Next, we complete the matrix with $n-p$ orthogonal vectors and find that the elements of our tangent space have the form $Z=X A+X_{\perp} B$ where $A$ is a skew-symmetric matrix $(A^T = -A)$.  


Introducing a metric
------
Until now, we considered the very general and abstract notion of a manifold. Since we want to optimize a function, we need to know the distance between two elements in order to move in the space and search for an optimum. We do this by defining an inner product on the <i>tangent space</i>. The first choice is the Euclidean inner product. For two elements of the tangent space $Z_1$ and $Z_2$ it is defined as $\left\langle Z_{1}, Z_{2}\right\rangle_{e}=\operatorname{tr}\left(Z_{1}^{T} Z_{2}\right)$. Using our representation from above and using the properties of $X$ an $X _{\perp}$, we expand this as follows: $$\langle Z, Z\rangle_{e} = tr (A^T A) + tr (B^T B) = \sum_{i>j} 2 a^{2}(i, j)+\sum_{i, j} b^{2}(i, j)$$
You see that the elements of $A$ are weighted twice as much as the elements of $B$. This is not desired, so instead we will use the <b>canonical inner product</b> $$\left\langle Z_{1}, Z_{2}\right\rangle_{c}=\operatorname{tr}\left(Z_{1}^{T}\left(I-\frac{1}{2} X X^{T}\right) Z_{2}\right)$$ It gives equal weights to both matrices.

Differentials and gradients
------
Next thing to consider when exploring the manifold are differentials. Since we want to use gradient-based optimization to solve the problem, we define the <i>differential</i> of a function $F : \mathcal{R}^{n \times p} \to \mathcal{R}$: $$D F_{X}(Z) =\sum_{i, j} \frac{\partial F}{\partial X_{i, j}} Z_{i, j} =\operatorname{tr}\left(G^{T} Z\right)$$ where $G=\left[\frac{\partial F}{\partial X_{i, j}}\right] \in \mathcal{R}^{n \times p}$

Under the canonical inner product, the vector $AX$ with $A = (GX^T − XG^T)$ represents the action of $DF_X$ on the tangent space $\mathcal{V}_{p}\left(\mathcal{R}^{n}\right)$.

TODO: refine and explain
Cayley Transform and the Search Curve
-----
The differential alone is not sufficient for finding an optimum. We need also a <i>descent curve</i>, which corresponds to the descent direction.


For this, let $W$ be any skew-symmetric matrix. Consider the curve $$Y(\tau)=\left(I+\frac{\tau}{2} W\right)^{-1}\left(I-\frac{\tau}{2} W\right) X$$ It has following properties:

1. It stays in the Stiefel manifold, i.e. $Y(t)^T Y(t) = I$
2. Its tangent vector at $\tau = 0$ is $Y'(0) = -WX$
3. If we set $W = A = GX^T − XG^T$, then the curve is a descent curve for $F$, i.e. the direction of steepest descent, which is also the negative gradient.


We can view $Y(\tau)$ as the point $X$ transformed by $\left(I+\frac{\tau}{2} W\right)^{-1}\left(I-\frac{\tau}{2} W\right)$. This is called the <i>Cayley transformation</i>

<b>Sketch of the algorithm:</b> Generate $X^{[k+1]}$ from $X^{[k]}$ by a curvilinear search along $Y(t)$ by changing $\tau$.

The search is carried out using the Armijo-Wolfe rule.

Since the matrix inversion is costly, we use the e Sherman-Morrison-Woodbury formula.

We get $Y(\tau)=X-\tau U\left(I+\frac{\tau}{2} V^{T} U\right)^{-1} V^{T} X$, where $U = [G, X]$ and $V = [X, -G]$ (this means concatenation).

This reduces the matrix size from $n \times n$ to $2p \times 2p$.

TODO: explain
Curvilinear search
-----
<b>Curvilinear Search:</b>
Initialize $\tau$ to a non-zero value.

$\begin{array}{l}
\text { Until }\left\{F(Y(\tau)) \leq F(Y(0))+\rho_{1} \tau F^{\prime}(Y(0))\right.\text { and } \\
\qquad \begin{array}{l}
\left.F^{\prime}(Y(\tau)) \geq \rho_{2} F^{\prime}(Y(0))\right\} \\
\text { do } \tau \leftarrow \frac{\tau}{2}
\end{array}
\end{array}$

Return $Y(\tau)$ as the curvilinear search "minimum".

To use curvilinear search, we need formulae for $F'(Y(τ))$ and $F'(Y(0))$.

$F^{\prime}(Y(\tau))=\operatorname{tr}\left(G^{T} Y^{\prime}(\tau)\right)$

where

$\begin{array}{l}
Y^{\prime}(\tau)=-\left(I+\frac{\tau}{2} A\right)^{-1} A\left(\frac{X+Y(\tau)}{2}\right) \\
Y^{\prime}(0)=-A X
\end{array}$

where, as before, $A=G X^{T}-X G^{T}$.

The resulting algorithm is a straightforward application of the curvilinear search.

<b>Minimize on Stiefel Manifold:</b>

Initialize: Set $k = 1$ and $X^{[k]}$ to a point in the manifold.

Iterate till convergence: 

1. Calculate $G^{[k]} = DF(X^{[k]})$, the $A$, $U$ and $V$ matrices.
2. Use curvilinear search to obtain $X^{[k+1]}$.

Test for convergence.
 

 Recap
 -----
 Let us sum up what we have done.
- Our approach:

    1. Find constraints of the problem
    2. Prove that the feasible set is a manifold
    3. Find the tangent space of the manifold
    4. Find/derive an inner product for the tangent space
    5. Find a representation of the gradient and with this the descent curve
    6. Use curvilinear search, i.e. gradient descent

You can use this strategy for similar non-convex constrained problems. The steps may differ or you may need additional steps depending on how easy it is to work in a specific manifold.

References
-----

[1]: W. M. Boothby, An Introduction to Differentiable Manifolds and Riemannian Geometry, Academic Press, 2002.
