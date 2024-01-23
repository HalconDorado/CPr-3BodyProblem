# Pseudocode #
```
= Initialize Variables =
Constants {
    - G: Gravitational constant
    - d: Earth Moon distance
    - Me: Earth Mass 
    - Mm: Moon Mass
    - P: Earth-Moon system orbital period
}
Variables/Lists {
    - E pos
    - M pos
    - ES distance
    - MS distance
    - Buffer of all above vars
    - aS: Static step size
    - aD: Dynamic step size
}
Vectors {
    - M position
    - M velocity
    - E position
    - E velocity
    - S position
    - S velocity
    - S acceleration
}

= Integrate =
While loop {
    - 
    - Change step size
    - Check termination condition
}
```

# Ideas #
What is the best size scale to use in order to minimize numerical precision error? \
\
More accurate representation of numbers by using more memory per variable -> potential tradeoff of speed vs. accuracy. \
\
Rotating reference frame or. Static reference frame?

# Milestone: #
Solve orbit of rocket at L2 point over 1 lunar orbit. Use both Taylor expansion & RK or SciPy.
Find optimal value of delta r (distance between L2-Moon) to error < 1% and plot orbit of
Moon & rocket. Determine separation of rocket-Moon after 1 orbit.

- adaptive step size methods: Romberg's Method and RKF Method
- 

# Extension: #
Galaxy collisions
- gas interactions -> fluid dynamics
- separation of the DM halo
- Approximation that M_dm >> M_stars?
- Initial population of disk according to observed distribution -> ``k1*exp(-z/z0)*k2*exp(-r/r0)``
- 

Kepler's Law things...

v=2pi*r/P
mv^2=GMm/r
4*pi^2*r^2/P^2=GM/r
P^2 = 4*pi^2*r^3/GM