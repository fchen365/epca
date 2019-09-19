# EPCA

EPCA: **e**xploratory **P**rinciple **C**omponent **A**nalysis (PCA) for data with inner nature of sparsity. It contains a series of implementations, including sparse \{matrix decomposition, PCA, independent component analysis (ICA)\} via **p**enalized-by-**r**otation **m**atrix **d**ecomposition (PRMD).


The PRMD is defined as the following optimization: 
\begin{eqnarray}
	&\min &  \| X-UDV^\top \|^2_F \label{eq:obj1}\\
	&\text{s.t.} & \mathcal{P_R}(U) \le \gamma_1, \mathcal{P_R}(V) \le \gamma_2 . \nonumber
%	&\text{s.t.} & \mathcal{P}(UO) \le \gamma_1, \mathcal{P}(VR) \le \gamma_2, \text{ for some } O,R\in\mathcal{U}(k).\nonumber 
\end{eqnarray}
where the penalizing operator is defined for any $$V\in\mathbb R ^{p \times k},V^\top V=I_k,$$
\begin{eqnarray*}
\mathcal{P_R}(V) = \min_{R\in\mathcal{U}(k)}\mathcal{P}_R(V), & \text{ where } & \mathcal{P}_R(V)=\mathcal{P}(VR).
\end{eqnarray*}