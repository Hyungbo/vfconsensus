We consider a multi-agent system given by

$$\dot x_i = [-x_i, \; \sin(t)] \; \theta_i + k \textstyle{\sum_{j \in \cN_i}} (x_j - x_i)$$

where $\theta_i \in \R^2$.
With each component of $\theta_i(0)$ varying within $[a,b]$ ($b>a>0$), both components of $\theta_0^* = (1/N) \sum_{i=1}^N \theta_i(0) =: [\alpha,\beta]^T$ are positive.
Then, it is straightforward to see that Assumptions 3 and 4 in the paper hold. 
Assumption 1 is also verified because equation (6) becomes $\dot {\hat s} = -\alpha \hat s + \beta \sin(t)$ and thus, $\hat s(t)$ tends to $A \sin(t + \gamma)$ with $A>0$ and $\gamma \not = 0$.

### Simulation

Number of agents N = 3.
Graph is the ring graph.
Coupling gain k = 1.



#### When there is no coupling
![When there is no coupling](toy1.png)

#### When there is coupling with no adaptation
![When there is coupling with no adaptation](toy1.png)

#### When there is coupling and adaptation
![When there is coupling and adaptation](toy1.png)
And, the progress of parameters
![Theta plot](theta.png)
