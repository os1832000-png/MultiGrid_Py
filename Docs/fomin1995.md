Pergamon

J. Quant. Spectrosc. Radiat. Transfer Vol. 53, No. 6, pp. 663–669, 1995

Copyright © 1995 Elsevier Science Ltd

0022-4073(95)00029-1

Printed in Great Britain. All rights reserved

0022-4073/95 $9.50 + 0.00

# EFFECTIVE INTERPOLATION TECHNIQUE FOR LINE-BY-LINE CALCULATIONS OF RADIATION ABSORPTION IN GASES

B. A. FOMIN

Russian Research Centre “Kurchatov Institute”, 123182 Moscow, Russia

(Received 20 July 1994)

Abstract—An interpolation technique for Line-by-Line calculations of absorption coefficients is developed. The technique gives a possibility considerably to increase the calculation speed. This work may be particularly useful in the calculation of radiation transfer in the atmosphere.

# INTRODUCTION

Selective gas absorption is a basic component of the interaction between the radiation and the atmosphere. Its accurate line-by-line (LBL) calculation is necessary for solution of a lot of problems of satellite meteorology, climate researches etc. because only the LBL approach provides the means of direct consideration of the detailed spectral structure of the gases from the microwave to the i.r.

The monochromatic volume absorption coefficient, $K(\nu)$, at any wavenumber $\nu$ can be calculated using the LBL technique by formula

$$
K(\nu) = \sum_{i} f_{i}(\nu, \tilde{\nu}_{i}) \tag{1}
$$

where $f_{i}(\nu, \tilde{\nu}_{i})$ is the absorption profile and $\tilde{\nu}_{i}$ is the position of the $i$th spectral line. The use of the LBL technique for calculation of a transfer of radiation generally leads to the integration over frequency

$$
\int_{\Delta \nu} K(\nu) \, \mathrm{d}\nu \tag{2}
$$

where $\Delta \nu$ may be $&gt;10^{4} \, \mathrm{cm}^{-1}$ in practice.

So the calculation from Eq. (1) generally needs to be performed at every points $\nu_{j}$ of a nonuniformly or uniformly spaced frequency grid as shown in Fig. 1.

The typical distance $H$ between the nodal points is $&lt; \sim 10^{-3} \, \mathrm{cm}^{-1}$ because the grid must be fine enough so that the narrowest line be adequately sampled$^{1}$ ($H \sim$ Doppler half width in the upper atmosphere).

The disadvantage of LBL calculations is that they require large computation times because a number of functions $f_{i}(\nu, \tilde{\nu}_{i})$ in sum (1) may reach tens of thousands. It would be noted that the use of the latest editions of the spectroscopic data bases for the LBL calculations requires more computation work than previous editions since they contain much more spectroscopic information (cf data bases HITRAN-86, 92$^{2,3}$). Moreover, for some problems the absorption coefficient $K(\nu)$ needs to be calculated over very large interval $\Delta \nu$, consequently, at a great number of points $\nu_{j}$ (to $\Delta \nu / H \sim 10^{4} / 10^{-3} \sim 10^{7}$).

So, the problem of effective LBL algorithms is of extreme importance and any improvements of the LBL technique may be of great practical usefulness although a lot of LBL computer codes have been developed since the sixties. Unfortunately, calculation details of some codes are not well described in the published papers. This makes the analysis of their efficiency difficult but not impossible. The point is that these codes generally use different interpolation techniques because
664
B. A. Fomin

![img-0.jpeg](img-0.jpeg)
Fig. 1. Illustration of the absorption spectrum $K(\nu_{j})$.

direct calculations by the formula (1) thousand times per a unit wavenumber interval (in $\mathrm{cm}^{-1}$) at points $\nu_{j}$ is prohibitive in practice. It should be noted that as a rule the calculations of functions $f_{i}(\nu_{j},\tilde{\nu}_{i})$ consume the main part of the computer work. On the contrary, the computer work of the interpolations may be negligible as explained below. So we can roughly evaluate the efficiency of any LBL code if we count the number of points where the absorption profiles must be calculated directly before the interpolations. In other words, we will compare the efficiencies of interpolations of the LBL algorithms only. (Of course the methods of the effective line shape calculation, of the line cut off, etc. are very important but they may be considered separately.)

Owing to the fact that the line profile changes slower at a distance from the line center than it does near the center (see Fig. 1) there are two basic methods which provide a possibility to increase the calculation speed. These methods have been utilized in most of the LBL codes since the well-known paper by Drayson⁴ was published. A starting point to the first of them is the second wavenumber grid to the far wings of lines. The typical distance between the nodal points of this grid is generally $1 - 10\mathrm{cm}^{-1}$ ($cfH \sim 0.001\mathrm{cm}^{-1}$) depending on the LBL code. At each point of the wider grid sum (1) must be calculated only over lines spaced at a greater distance from it than $D \sim 3 - 7\mathrm{cm}^{-1}$ (see Refs. 4–6 etc.). Then it is only necessary to perform a simple polynomial interpolation for intermediate wave numbers $\nu_{j}$ of the first fine grid to consider all these far-off lines. The computer time for the interpolation (generally quadratic) is much about the same as for the calculation of only one Lorentz line profile and hence is insignificant. As a result, using the first method the far-end wing contributions may be very effectively allowed for by the precalculated ‘continuum’ at the points of the wider grid. So the main problem now is how to take into account the contributions from the closer lines. In other words, the problem is to carry out summation (1) at a lot of $\nu_{j}$ points over $i$th lines if $|\nu_{j} - \tilde{\nu}_{i}| &lt; D$. Thus, as it can be easy seen, every $i$th line profile $f_{i}(\nu_{j},\tilde{\nu}_{i})$ needs to be calculated at $2D / H \sim 10 / 10^{-3} \sim 10^{4}\nu_{j}$ mesh points.

The second basic method which has been used in several LBL codes just gives a possibility to decrease a number of these points and in this way to decrease the computer time. The principal feature of the method is nonuniformly spaced wavenumber grid instead of uniform grid $\nu_{j}$. The mesh points of a nonuniform grid are clustered near each line centre $\tilde{\nu}_{i}$ to attain the acceptable accuracy of integration (2) by a corresponding quadrature (generally the Legendre–Gauss quadrature) which is based on the interpolation. Therefore, the mesh interval increases between the lines starting from $H$ at their centers (the details see in Refs. 4–6, etc.). It is clear that this method is effective when the distances between the lines are considerably greater than the typical line half width $\alpha$ ($\alpha \approx H$). But for the most time-consuming cases when the density of the lines is high (in strong bands) a number of mesh points of both uniform and nonuniform grids becomes approximately the same and, consequently, the second basic method is not very effective. Moreover, a nonuniform grid is not as easily applicable in LBL codes as a uniform grid. So, the majority of contemporary LBL codes use a uniform grid (GENLN2 Ref. 6, etc.).

In the final analysis all LBL codes which use single grid only (except for the grid of the continuum) have a common fundamental disadvantage because they calculate each $i$th absorption profile at the mesh points clustered near the centers $\dots, \tilde{\nu}_{i-2}, \tilde{\nu}_{i-1}, \tilde{\nu}_{i+1}, \tilde{\nu}_{i+2}, \dots$ of the other lines
LBL calculations of radiative absorption

although for the $i$th line a fine grid is necessary at center $\tilde{\nu}_i$ only. The technique which gives a possibility to avoid this disadvantage has been applied in the widely distributed LBL code FASCOD. $^{7}$ The principal feature of this code is a decomposition of each line profile into three functions "fast", "slow" and "very slow" which then must be considered individually using three uniform grids with the wavenumber steps $\alpha/4$, $\alpha$ and $4\alpha$, respectively. Using this technique only 99 mesh points are required to describe Voight profile over $(\tilde{\nu} \pm 64\alpha)$ interval (see Ref. 7) instead of $2 * 64\alpha / 0.25\alpha = 512$ points, but it takes to overcome some difficulties concerning the decomposition of the line profile. The difficulties may be considerable if it is necessary to take into account the line mixing effects when we have no Voight profile. The simple interpolation technique which gives possibility to avoid the decomposition is given in the next Sec. Moreover, using the technique suggested one may calculate each absorption profile over $(\tilde{\nu} \pm 64\alpha)$ interval at $\sim 30-50$ points only instead of 99 points and this is the way still more to decrease the computer time by 2-3 times.

## THE LINE-BY-LINE TECHNIQUE

At first we shall consider the interpolation problem for a single line only in the region $D &gt; |\nu - \tilde{\nu}| &gt; T$ around the line center $\tilde{\nu}$ where $T = (2 \div 3)(\alpha_L + \alpha_D)$ and $\alpha_L / \alpha_D$ is a Lorentz/Doppler half width. We are not going to scrutinize this problem near the center $|\nu - \tilde{\nu}| &lt; T$ since it may be easily solved by a simple polynomial interpolation using approximately a tenth of mesh points of the uniform grid. The use of the nonuniform grid instead of the uniform one is not very useful in this case because the number of mesh points there is commonly quite restricted. But the nonuniform grid is very useful to obtain any effective interpolating technique for the complete wing of the line. For this aim both left and right wings (except for the central region) must be partitioned into the $2L$ portions located at unequal intervals $T &lt; |\nu - \tilde{\nu}| &lt; t_1, t_1 &lt; |\nu - \tilde{\nu}| &lt; t_2, \ldots, t_{L-1} &lt; |\nu - \tilde{\nu}| &lt; t_L = D$ where $T, t_1, t_2, \ldots, t_L$ are distances between the line center and the boundaries of the intervals. Factor 2 arises owing to the right and left wings. We shall assume that a common interpolating formula and the same number of the mesh points will be used for each interval as it is in the case of majority of LBL codes.

In order to obtain the effective interpolating technique we shall use the fact that any spectral line at $|\nu - \tilde{\nu}| &lt; 3\mathrm{cm}^{-1} \sim D$ has a shape which is quite close to the direct proportion to $(\nu - \tilde{\nu})^{-2}$ if $|\nu - \tilde{\nu}| &gt; T$. (The correction factor can be evaluated by formula $[1 + (1.5\alpha_D^2 - \alpha_L^2) / (\nu - \tilde{\nu})^2]$. This formula can be easily obtained if one considers the Voight integral when $\alpha_L / |\nu - \tilde{\nu}| \ll 1$ and $\alpha_D / |\nu - \tilde{\nu}| \ll 1^8$).

Owing to this fact we define the intervals by $t_1 = C * T$, $t_2 = C^2 * T, \ldots, t_L = C^L * T$ where $C &gt; 1$ is a parameter which will be considering below. It gives relative interpolating errors of function $(\nu - \tilde{\nu})^{-2}$ (i.e. the line wing) identical in each of the intervals. This statement may be easily checked if we consider the Lagrange's interpolating formula for this function in the above intervals. Parameter $C$ may be different in different codes. For example, Aoki has used $C = 1.3$. But in the present work we shall use $C = 2.0$ to obtain intervals which are multiply of the central $(\tilde{\nu}, \tilde{\nu} \pm T)$ interval. This property will be used below. Thus we have intervals

$$
(\tilde {\nu} - D, \tilde {\nu} - 2 ^ {L - 1} T), \dots , (\tilde {\nu} - 2 ^ {2} T, \tilde {\nu} - 2 T), (\tilde {\nu} - 2 T, \tilde {\nu} - T),
$$

$$
(\tilde {\nu} + T, \tilde {\nu} + 2 T), (\tilde {\nu} + 2 T, \tilde {\nu} + 2 ^ {2} T), \dots , (\tilde {\nu} + 2 ^ {L - 1} T, \tilde {\nu} + D) \tag {3}
$$

where the interpolating formula should be applied. The interpolation order and, consequently, the number of interpolating nodal points may be changed to reach the acceptable accuracy. For example, the quadratic interpolation formula using three equidistant points $(\tilde{\nu} \pm T, \tilde{\nu} \pm 1.5T, \tilde{\nu} \pm 2T)$ in the closer intervals to the center $\tilde{\nu}$ and $(\tilde{\nu} \pm 2T, \tilde{\nu} \pm 3T, \tilde{\nu} \pm 4T)$ for the next intervals, etc. gives us the worst errors $\delta = 7.8\%$ at points $\tilde{\nu} \pm 1.8T, \tilde{\nu} \pm 3.6T, \ldots$. But the interpolation using five equidistant points $(\tilde{\nu} \pm T, \tilde{\nu} \pm 1.25T, \tilde{\nu} \pm 1.5T, \tilde{\nu} \pm 1.75T, \tilde{\nu} \pm 2T$ etc.) gives us the worst errors $\delta = 0.4\%$ at points $\tilde{\nu} \pm 1.91T, \tilde{\nu} \pm 3.82T, \ldots$. All these results have been obtained by simple numerical experiments for function $|\nu - \tilde{\nu}|^{-2}$ (it is enough to consider only one interval). It should be stressed out that this function is much poorer 'interpolatable' than any Lorentz shape and generally Voight shape owing to its much faster variations. So, as a matter of experience, the errors must be less than $\delta$ in the closer intervals to $\tilde{\nu}$.
B. A. Fomin

The number of the intervals for the each wing may be easily evaluated by the equation

$$
2^{L} * T = D. \tag{4}
$$

It takes only a dozen intervals to approximate the wing starting with $T = 0.001\,\mathrm{cm}^{-1}$ to $D = 4\,\mathrm{cm}^{-1}$ ($2^{12} * 0.001 = 4.096$). Thus, the absorption profile needs to be calculated directly only at $2(12 + 12) \sim 50$ mesh points to approximate both wings with accuracy better than $7.8\%$ by the quadratic interpolation or $(12 + 12) * 4 \sim 100$ points with accuracy $0.4\%$ by the five-point interpolation (one point is shared by neighboring intervals). There is an opportunity to decrease the number of such intervals by changing parameter $C$ from $C = 2$ to $C = 3$, $C = 4$, etc. which also gives us the intervals multiple of the central ($\vec{v}, \vec{v} \pm T$) one. For example, if $C = 4$ the number of intervals equals a half of dozen only. But numerical experiments prove that using $C = 3$, $C = 4$, ... is not worthwhile because it leads to increase of the errors or needs a greater number of the mesh points in each interval. For instance, for both cases, when $C = 4$ with the five-point interpolation and $C = 2$ with quadratic interpolation, the same number of mesh points is required but they give $\delta \sim 10\%$ and $\delta \sim 8\%$, respectively. In the first case the optimization of distances between the mesh points has been done. But the use for the first case the same points that for the second case (here must be considered two neighbor intervals) gives $\delta = 66\%$. As a result, we shall use the doubling wavenumber intervals [Eq. (3)] only. The interpolation order may be 2, 3, 4, ... etc. as stated above. But in the present work we shall use mainly the quadratic formula to simplify further explanation of the method. Moreover, it gives the accuracy which is acceptable for majority of practical applications. Thus, one can to calculate any absorption profile by $\sim 100$ times more rapidly than when the uniform grid is used (cf. $\sim 2D / H \sim 10^3 + 10^4$ and 50 points). It will be interesting to compare this technique with the others which have been described in the papers Refs. 4–6, 9 and etc. Unfortunately, there are no calculation details to compare the errors and the total number of the mesh points. But the comparison with the FASCOD$^7$ (one of the most effective codes) may be easily performed. As above mentioned, the FASCOD calculates the line profile at 99 points over the interval ($\vec{v} \pm 64\alpha$) using three uniform grids with the steps $\alpha/4$, $\alpha$ and $4\alpha$. In our case it takes only 33 points

$$
\left(\vec{v}, \vec{v} \pm \frac{\alpha}{4}, \vec{v} \pm \frac{\alpha}{2}, \vec{v} \pm \frac{3}{4}\alpha, \vec{v} \pm \alpha, \vec{v} \pm \frac{3}{2}\alpha, \vec{v} \pm 2\alpha, \vec{v} \pm 3\alpha, \dots, \vec{v} \pm 48\alpha, \vec{v} \pm 64\alpha\right)
$$

for the quadratic interpolation ($\delta \sim 7.8\%$) or 57 points

$$
(\dots, \vec{v} \pm \alpha, \vec{v} \pm \frac{3}{4}\alpha, \vec{v} \pm \frac{3}{2}\alpha, \vec{v} \pm \frac{7}{4}\alpha, \vec{v} \pm 2\alpha, \dots, \vec{v} \pm 56\alpha, \vec{v} + 64\alpha)
$$

for the five-point interpolation ($\delta \sim 0.4\%$). Thus, this interpolation technique is more effective by 2–3 times.

Now we are ready to consider the interpolation technique for the case of a great number of lines. The principal feature of the technique is a series of the uniform grids $v_{j}^{(l)}$ with the doubling wavenumber steps $h_{i}$:

$$
h_{0} = H, \quad h_{1} = H * 2, \quad h_{2} = h_{1} * 2, \quad \dots, \quad h_{l} = H * 2^{l}, \quad \dots \quad l = 0, 1, \dots, L
$$

$$
v_{j}^{(l)} = v_{\text{start}} + h_{i} * j \quad j = 0, 1, \dots \tag{5}
$$

where $v_{j}^{(0)} = v_{j}$ is a fine grid, $v_{\text{start}}$ is a starting wavenumber, $L$ is the number of the wider grids ($L \approx 10$).

Every $l$th grid is desired for the above-mentioned $l$th interpolating interval [Eq. (3)] of each $i$th line as it is shown in Fig. 2. By means of this series, one may consider any $i$th line independently on another one using $v_{j}^{(l)}$ points to obtain the above mentioned individual effective nonuniform grid. This trick solves the problem. It is similar to the above-mentioned Drayson’s idea$^4$ which has used the second grid for far line contributions only.

Some remarks should be made. Due to the fact that $v_{j}^{(l)}$ mesh points are fixed the distances between the line centers and the boundaries of the interpolation intervals may be greater than required for relations [Eq. (3)] but cannot be less than that in order to avoid the uncontrolled errors. So two neighboring interpolation intervals may be located at the same grid as it is shown in Fig. 2.
LBL calculations of radiative absorption
667

![img-1.jpeg](img-1.jpeg)
Fig. 2. Illustration of the interpolation technique using a series of the grids [Eq. (5)].

As a whole the grid using $v_{j}^{(l)}$ points has the properties similar to those of grid [Eq. (3)] where the parameter $C = 2$ has been used.

It should be stressed out that summation (1) must be taken individually for each grid $v_{j}^{(l)}$ with respect to the corresponding portion of the $i$th line shape before interpolation

$$
\tilde{\varphi}_{j}^{(l)} = \tilde{\varphi}_{j}^{(l)} + f_{i}(v_{j}^{(l)}, \tilde{v}_{i}),
$$

$$
\varphi_{j+1}^{(l)} = \varphi_{j+1}^{(l)} + f_{i}(v_{j+1}^{(l)}, \tilde{v}_{i}),
$$

$$
\tilde{\tilde{\varphi}}_{j+2}^{(l)} = \tilde{\tilde{\varphi}}_{j+2}^{(l)} + f_{i}(v_{j+2}^{(l)}, \tilde{v}_{i}), \quad j = 0, 2, \dots \tag{6}
$$

where $\tilde{\varphi}_j^{(l)},\varphi_{j + 1}^{(l)},\tilde{\tilde{\varphi}}_{j + 2}^{(l)}$ are accumulated contributions at $v_{j}^{(l)},v_{j + 1}^{(l)},v_{j + 2}^{(l)}$ points which will be used as three nodal points. The values $\tilde{\varphi}_j^{(l)}$ and $\tilde{\tilde{\varphi}}_j^{(l)}$ at each even point $v_{j}^{(l)}$ must be considered individually (it takes two storage locations) depending on whether the point is on the left or on the right of the central odd point $v_{j + 1}$ in [Eq. (6)] to take into account the line shape "breaks" which are shown in Fig. 2. It is very convenient that

$$
v_{0}^{(l-1)} = v_{0}^{(l)}, \quad v_{1}^{(l-1)} = 0.5(v_{0}^{(l)} + v_{1}^{(l)}), \quad v_{2}^{(l-1)} = v_{1}^{(l)},
$$

$$
v_{3}^{(l-1)} = 0.5(v_{1}^{(l)} + v_{2}^{(l)}), \dots \tag{7}
$$

where $l = 2,3,\ldots ,L$

Moreover, if $T = 4H$ it is easy to see that $v_{0} = v_{0}^{(l)}, v_{1} = 0.5(v_{0}^{(1)} + v_{1}^{(1)}), v_{2} = v_{1}^{(1)}, \ldots$. So we shall use a quadratic interpolation of any function $p(\nu)$ at middle points between the nodals only

$$
p(\nu^*) \approx 0.375p(v_{a}) + 0.75p(v_{b}) - 0.125p(v_{c})
$$

$$
p(\nu^{**}) \approx 0.375p(v_{c}) + 0.75p(v_{b}) - 0.125p(v_{a}) \tag{8}
$$

where

$$
v_{b} = 0.5(v_{a} + v_{c})
$$

$$
v^* = 0.5(v_{a} + v_{b})
$$

$$
v^{**} = 0.5(v_{c} + v_{b}).
$$

When all the lines are taken into account, then using Eqs. (7) and (8) one may obtain absorption coefficient $K(v_{j})$ at every point $v_{j}$ of the fine grid using the recurrence procedure

$$
\tilde{\varphi}_{m}^{(l-1)} = \tilde{\varphi}_{m}^{(l-1)} + \tilde{\varphi}_{m/2}^{(l)}
$$

$$
\varphi_{m+1}^{(l-1)} = \varphi_{m+1}^{(l-1)} + 0.375 * \tilde{\varphi}_{m/2}^{(l)} + 0.75 * \varphi_{m/2+1}^{(l)} - 0.125 * \tilde{\tilde{\varphi}}_{m/2+2}^{(l)},
$$

$$
\tilde{\varphi}_{m+2}^{(l-1)} = \tilde{\varphi}_{m+2}^{(l-1)} + \varphi_{m/2+1}^{(l)}, \quad \tilde{\tilde{\varphi}}_{m+2}^{(l-1)} = \tilde{\tilde{\varphi}}_{m+2}^{(l-1)} + \varphi_{m/2+1}^{(l)}
$$
668
B. A. Fomin

![img-2.jpeg](img-2.jpeg)
Fig. 3. Lorentz line shape normalized to unity at the line center $[f_{i}(\tilde{\nu}_{i},\tilde{\nu}_{i}) = 1]$ calculated directly (solid line) and by the interpolation [Eq. (9)] (dashed line), $\alpha /h = 1$.

$$
\varphi_{m + 3}^{(l - 1)} = \varphi_{m + 3}^{(l - 1)} - 0.125 * \tilde{\varphi}_{m / 2}^{(l)} + 0.75 * \varphi_{m / 2 + 1}^{(l)} + 0.375 * \tilde{\tilde{\varphi}}_{m / 2 + 2}^{(l)},
$$

$$
\tilde{\tilde{\varphi}}_{m + 4}^{(l - 1)} = \tilde{\tilde{\varphi}}_{m + 4}^{(l - 1)} + \tilde{\tilde{\varphi}}_{m / 2 + 2}^{(l)} \tag{9}
$$

where $m = 0,4,8,\ldots ,l = L,L - 1,\ldots ,2$

$$
K \left(v _ {m + 1}\right) = K \left(v _ {m + 1}\right) + 0.375 * \tilde {\varphi} _ {m / 2} ^ {(l)} + 0.75 * \varphi_ {m / 2} ^ {(l)} - 0.125 * \tilde {\tilde {\varphi}} _ {m / 2 + 2} ^ {(l)}
$$

$$
K \left(v _ {m + 2}\right) = K \left(v _ {m + 2}\right) + \varphi_ {m / 2 + 1} ^ {(1)}
$$

$$
K \left(v _ {m + 3}\right) = K \left(v _ {m + 3}\right) - 0.125 * \tilde {\varphi} _ {m / 2} ^ {(l)} + 0.75 * \varphi_ {m / 2 + 1} ^ {(l)} + 0.375 * \tilde {\tilde {\varphi}} _ {m / 2 + 2} ^ {(l)}
$$

$$
K \left(v _ {m + 4}\right) = K \left(v _ {m + 4}\right) + \tilde {\varphi} _ {m / 2 + 2} ^ {(1)}.
$$

It is very important that this simple procedure is performed only once for it to consume the negligible computer time.

A few more words should be said about the precision. The typical numerical errors of the calculation of the line shape by this technique are shown in Fig. 3. It is clearly seen that the errors from different contributions $f_{i}(\nu, \tilde{\nu}_{i})$ (they are less than $\delta \sim 7.8\%$ as it has been obtained above) are neutralized partially in sum (1) [i.e. in $K(\nu)$] because they alternate their signs. Furthermore, it is useful to consider error $\delta^{*}$ which is the maximum error divided by the line shape amplitude $f_{i}(\tilde{\nu}_{i}, \tilde{\nu}_{i})$ and relative errors $\Delta$ of the numerical integration (using the parabolic rule) over the interval. The first $\delta^{*}$ is $&lt; 1\%$ if $\alpha_{L} / H &gt; 1$ as one may obtain from a simple numerical experiment (see Ref. 10 as well). The use of function $|\nu - \tilde{\nu}_{i}|^{-2}$ gives us the errors on the safe side $\Delta = 1 / 108 &lt; 1\%$. In the final analysis this technique, as a rule, gives us the precision better than $1\%$. This precision may be easily improved if the five-point interpolation is used instead of a quadratic one. For this aim it is necessary only to insert point $0.5 * (\nu_{j}^{(l)} + \nu_{j+1}^{(l)})$ between $\nu_{j}^{(l)}, \nu_{j+1}^{(l)}$ points in every $l$th grid.

Generally speaking, the uniform fine grid $\nu_{j}$ may not be effective for further computations of radiative transfer parameters. But it does not make any problem to obtain $K(\nu)$ at the points of a more effective nonuniform grid if prior to interpolation [Eq. (9)] one omits $\nu_{j}$ points where $K(\nu_{j}) = 0$, then $\nu_{j}^{(l)}$ points where $\varphi_{j}^{(l)} = 0$ etc. and performs the interpolations [Eq. (9)] for the remaining points only.

Some remarks should be made about the computer memory and time consumption. If the fine grid has $n_0$ points, one may obtain the total number of points $N$ using [Eq. (3)]

$$
\begin{array}{l}
N = n _ {0} + n _ {1} + n _ {2} + \dots + n _ {L} \approx n _ {0} + \frac {3}{4} (n _ {0} - 1) + \frac {3}{8} (n _ {0} - 1) + \dots + \frac {3}{2 * 2 ^ {L}} (n _ {0} - 1) \\
= n _ {0} + \frac {3}{2} \left(1 - \frac {1}{2 ^ {L}}\right) * (n _ {0} - 1) &lt; 2.5 n _ {0}. \tag{10}
\end{array}
$$

where $n_1, \ldots, n_L$ are the numbers of the storage locations for every grid.

Using the above technique each line shape must be calculated directly at $\sim 50$ points only. It is evident that the number of these points increases very slowly (as logarithm to base 2 instead of the inversely proportional dependence) if the step of fine grid $H$ decreases. The computer time for
LBL calculations of radiative absorption

the recurrence procedure (9) is negligible. Finally the proposed algorithm is in need of the double-triple computer storage but gives a possibility significantly to speed up the computation work.

# CONCLUSIONS

Using this technique the effective software has been developed to obtain the benchmark calculations of the solar radiation in the cloudy and aerosol atmospheres for climate problems. Moreover, the software has been used for some satellite sounding problems. An average speed of $\sim 10^{3}\mathrm{cm}^{-1}/\mathrm{h}$ has been obtained for usual LBL calculations of the absorption coefficients ($\sim 50$ levels and the spectral resolution of $\sim 0.5*10^{-3}\mathrm{cm}^{-1}$) using an IBM PC 486 computer.

Acknowledgements—I would like to acknowledge the invaluable support received from the Institute of Molecular Physics of Russian Research Centre “Kurchatov Institute” where the work was performed. I am grateful to L. S. Rothman (GL/OPI, USA) and other members of ASA (Atmospheric Spectroscopy Applications) working group for providing me with spectroscopic data bases. In particular I would like to thank S. A. Clough (Atmospheric and Environmental Research Inc., U.S.A.) and D. P. Edwards (NCAR, USA) for sending me the extremely useful and well-documented information on their line-by-line models. Finally, I wish to thank A. N. Rublev, A. N. Trotsenko, Sh. Sh. Nabiev and other colleagues from the “Kurchatov Institute” for many useful discussions.

# REFERENCES

1. D. P. Edwards and L. L. Strow, J. Geophys. Res. 96, 20859 (1991).
2. L. S. Rothman, R. R. Gamache, A. Goldman, L. R. Brown, R. A. Toth, H. M. Pickett, R. L. Poynter, J.-M. Flaud, C. Camy-Peyret, A. Barbe, N. Husson, C. P. Rinsland, and A. H. Smith, Appl. Opt. 26, 4058 (1987).
3. L. S. Rothman, R. R. Gamache, R. Tipping, C. P. Rinsland, M. A. H. Smith, D. C. Benner, V. Malathy Devi, J.-M. Flaud, C. Camy-Peyret, A. Goldman, S. T. Massie, L. R. Brown, and R. A. Toth, JQSRT 48 (1992).
4. S. R. Drayson, Appl. Opt. 5, 385 (1966).
5. V. G. Kunde and W. C. Maguire, JQSRT 14, 803 (1974).
6. D. P. Edwards, NCAR Technical Note. NCAR/TN-367 + STR, Boulder, CO, p. 35 (1992).
7. H. J. P. Smith, D. J. Dube, M. E. Garden, S. A. Clough, F. X. Kneizys, and L. S. Rothman, Rep. Air Force Geophys. Lab. AFGL-TR-78-0081, Hanscom, MA (1978).
8. J. Humlicek, JQSRT 27, 437 (1982).
9. T. Aoki, Papers Meteor. Geophys 39, 53 (1988).
10. B. A. Fomin, Preprint of “Kurchatov Institute”, IAE-5658/1, Moscow, (1993).