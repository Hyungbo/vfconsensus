Consider a multi-agent system given by
$$\dot x_i = [-x_i, \; \sin(t)] \; \theta_i + k \textstyle{\sum_{j \in \cN_i}} (x_j - x_i)$$
where $\theta_i \in \R^2$.
With each component of $\theta_i(0)$ varying within $[a,b]$ ($b>a>0$), both components of $\theta_0^* = (1/N) \sum_{i=1}^N \theta_i(0) =: [\alpha,\beta]^T$ are positive.
Then, it is straightforward to see that Assumptions \ref{ass:psi} and \ref{ass:bd} hold. 
Assumption \ref{ass:pe} is also verified because \eqref{eq:bd} becomes $\dot {\hat s} = -\alpha \hat s + \beta \sin(t)$ and thus, $\hat s(t)$ tends to $A \sin(t + \gamma)$ with $A>0$ and $\gamma \not = 0$.
Fig.~\ref{fig:2} shows convergence of parameters (when $N=3$).
