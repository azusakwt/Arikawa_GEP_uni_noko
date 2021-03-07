# 有川湾　GEP　推定

# The model


**PPFD Model**

\begin{aligned}

log(\mu_p) &= f(m) + f(m | state) + state \\
\alpha &= \mu / \beta \\
y_p &\sim Ga(\alpha, \beta) \\
\beta &\sim N^{+}(0, 1)\\

\end{aligned}

**Temperature Model**

$$
\mu_t = f(m) + f(m | state) + f(m | location) + state + location
$$

**GEP Model**

$$
\mu_g = y_p + y_t + state + location
$$


