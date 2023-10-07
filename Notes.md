# Pseudocode #
```
= Initialize Variables =
Constants {
    - G: Gravitational constant
    - d: Earth Moon distance
    - M_E: 
    - M_M: 
}
Variables {
    - E pos
    - M pos
    - ES distance
    - MS distance
    - a: Step size (dynamic or static?)
}
Vectors {
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
More accurate representation of numbers by using more memory per variable -> potential tradeoff of speed vs. accuracy.

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