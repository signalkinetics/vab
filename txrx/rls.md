# Matlab-style RLS update function.

Formatting: $v$ is $1\times1$, $\mathbf{v}$ is $n \times 1$, and $\mathbf{V}$ is $n \times n$, note that$\mathbf{v}^H \mathbf{V} \mathbf{v}$ is also $1\times1$


Loop invariant:
$$
\begin{aligned}
	\mathbf{R}_{uu} \mathbf{P} &= \mathbf{I}\\
	\mathbf{c} &= \mathbf{P} \mathbf{r}_{ud}
\end{aligned}
$$

Basic update rule:
$$
\begin{aligned}
\mathbf{R}_{uu}^\mathit{new} &=  \lambda\mathbf{R}_{uu} + \mathbf{u} \mathbf{u}^H \\

\mathbf{r}_{ud}^\mathit{new} &=  \lambda\mathbf{r}_{ud} + d^* \mathbf{u} \\
\end{aligned}
$$


## Validation of actual update rule.

Kalman gain is defined:
$$
\mathbf{k} = \frac{\mathbf{P}\mathbf{u}}{\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}} 
$$

Update of P:
$$
\begin{aligned}
	\mathbf{P}^\mathit{new} &= \frac{\mathbf{P} -  \mathbf{k}\mathbf{u}^H \mathbf{P}}{\lambda}\\

													&\mathbf{R}_{uu}^\mathit{new}\mathbf{P}^\mathit{new}   \\
													&= 

	\left(\lambda\mathbf{R}_{uu} + \mathbf{u} \mathbf{u}^H\right) \left(\frac{\mathbf{P} -  \mathbf{k}\mathbf{u}^H \mathbf{P}}{\lambda}\right) \\
													&= 
													\frac{1}{\lambda\left(\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}\right)} 
	\left(\lambda\mathbf{R}_{uu} + \mathbf{u} \mathbf{u}^H\right) \left(({\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}})\mathbf{P} -  \mathbf{P}\mathbf{u}\mathbf{u}^H \mathbf{P}\right) \\
													&=

													\frac{1}{\lambda\left(\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}\right)} 
	\left(
	({\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}})\lambda\cancel{\mathbf{R}_{uu} \mathbf{P} } +
	({\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}}) \mathbf{u} \mathbf{u}^H \mathbf{P}
-  \lambda\cancel{\mathbf{R}_{uu} \mathbf{P}}\mathbf{u}\mathbf{u}^H \mathbf{P}
-  \mathbf{u} \mathbf{u}^H \mathbf{P}\mathbf{u}\mathbf{u}^H \mathbf{P}
\right) \\
													&=
													\frac{1}{\lambda\left(\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}\right)} 
	\left(
	({\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}})\lambda \mathbf{I} +
	({\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}}) \mathbf{u} \mathbf{u}^H \mathbf{P}
-  \lambda \mathbf{u}\mathbf{u}^H \mathbf{P}
-  \mathbf{u} \mathbf{u}^H \mathbf{P}\mathbf{u}\mathbf{u}^H \mathbf{P}
\right) \\

													&= \mathbf{I} + 
													\frac{1}{\lambda\left(\lambda + \mathbf{u}^H\mathbf{P} \mathbf{u}\right)} 
	\left(
	\cancel{\lambda  \mathbf{u} \mathbf{u}^H \mathbf{P}
-  \lambda \mathbf{u}\mathbf{u}^H \mathbf{P}}
\cancel{	+ (\mathbf{u}^H\mathbf{P} \mathbf{u}) \mathbf{u} \mathbf{u}^H \mathbf{P}
-  \mathbf{u} (\mathbf{u}^H \mathbf{P}\mathbf{u})\mathbf{u}^H \mathbf{P}}
\right)  \\
													&= \mathbf{I}
\end{aligned}
$$

Update of c:
Define the equalization error:
$$
e = \mathbf{c}^H \mathbf{u} - d = \mathbf{r}^H \mathbf{P} \mathbf{u} - d
$$

$$
\begin{aligned}
	\mathbf{c}^\mathit{new} - \mathbf{c} 
	&=  \mathbf{P}^\mathit{new}\mathbf{r}_{ud}^\mathit{new} - \mathbf{P}\mathbf{r}_{ud}\\
	&= \left(\frac{\mathbf{P} -  \mathbf{k}\mathbf{u}^H \mathbf{P}}{\lambda}\right)
	\left(\lambda\mathbf{r}_{ud} + d^* \mathbf{u}\right) 
	- \mathbf{P}\mathbf{r}_{ud} \\ 
	&=

	\cancel{\frac{\mathbf{P}}{\lambda} \lambda\mathbf{r}_{ud}}
	 -  \frac{\mathbf{k}\mathbf{u}^H \mathbf{P}}{\lambda}
d^* \mathbf{u}
	-\mathbf{k}\mathbf{u}^H \mathbf{P} \mathbf{r}_{ud} 

	+  \frac{\mathbf{P}}{\lambda}
	d^* \mathbf{u}
	\cancel{- \mathbf{P}\mathbf{r}_{ud}} \\ 
	&= - \frac1{\lambda}(
	d^*  (\mathbf{u}^H \mathbf{P} \mathbf{u}) \mathbf{k} + \lambda  \mathbf{k}\mathbf{u}^H \mathbf{P} \mathbf{r}_{ud} -  
	d^* \mathbf{P} \mathbf{u}
	)\\


	&= - \frac1{\lambda}(
	d^*  (\mathbf{u}^H \mathbf{P} \mathbf{u}) \mathbf{k} 
	+ \lambda  \mathbf{k} d^*
	+ \lambda  \mathbf{k}\mathbf{u}^H \mathbf{P} \mathbf{r}_{ud} 
	- \lambda  \mathbf{k} d^*
	-d^* \mathbf{P} \mathbf{u}
	) \\

	&= - \frac1{\lambda}(
	\cancel{ (\mathbf{u}^H \mathbf{P} \mathbf{u} + \lambda ) \mathbf{k} d^*}
	+ \lambda  \mathbf{k}( \mathbf{u}^H \mathbf{P} \mathbf{r}_{ud}  -  d^*)
	\cancel{-d^* \mathbf{P} \mathbf{u}}
	)\\
	&=-\mathbf{k} e^*
\end{aligned}
$$
