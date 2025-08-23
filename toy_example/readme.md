# Example: First order system 

We consider a multi-agent system given by

$$\dot x_i = [-x_i, \sin(t)] \; \theta_i + k \textstyle{\sum_{j \in N_i}} (x_j - x_i)$$

where $\theta_i \in R^2$.
With each component of $\theta_i(0)$ varying within $[a,b]$ ($b>a>0$), both components of 

$$\theta_0^* = (1/N) \sum_{i=1}^N \theta_i(0) =: [\alpha,\beta]^T$$

are positive.
Then, it is straightforward to see that Assumptions 2 and 3 in the paper hold. 
Assumption 4 is also verified because equation (9) becomes 

$$\dot {\hat s} = -\alpha \hat s + \beta \sin(t)$$ 

and thus, $\hat s(t)$ tends to $A \sin(t + \gamma)$ with $A>0$ and $\gamma \not = 0$.

### Simulation

Number of agents N = 3.
Graph is the ring graph.
Coupling gain k = 1.

Initial states $x(0) = [x_1(0), x_2(0), x_3(0)]^T$ is taken as
x0 = rand(N,1);

And the initial $\theta_i(0)$ is taken as

$\theta_1(0) = [1, 0.9]^T$, $\theta_2(0) = [1.1, 0.8]^T$, $\theta_3(0) = [1.2, 0.7]^T$

#### When there is no coupling
![When there is no coupling](toy1.png)

#### When there is coupling with no adaptation
![When there is coupling with no adaptation](toy2.png)

#### When there is coupling and adaptation
![When there is coupling and adaptation](toy3.png)

And, the progress of parameters
![Theta plot](theta.png)


# Example: Second order system 

The paper considers the first order system only, but the proposed adaptation law still works for higher order systems.
This time we consider

$$
\begin{aligned}
\dot x_i^1 &= [-x_i^1, \sin(t), 1] \cdot \theta + k \textstyle{\sum_{j \in N_i}} (x_j^1 - x_i^1) \\
\dot x_i^2 &= [\sin(2t), -x_i^2, 1] \cdot \theta + k \textstyle{\sum_{j \in N_i}} (x_j^2 - x_i^2)
\end{aligned}
$$

where $\theta \in R^3$.

The simulation code for this case is `toy_sim2.m`.

