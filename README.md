# Gradient-Based Optimizer (GBO) and RUNge Kutta optimizer (RUN) for Reservoir Operation Optimization

## Introduction

This repository implements the **Gradient-Based Optimizer (GBO)**, a general purpose novel metaheuristic optimization algorithm inspired by gradient search techniques. The GBO algorithm was introduced by Ahmadianfar et al. in their 2020 paper: ([Gradient-Based Optimizer: A New Metaheuristic ... ]([Ahmadianfar, Iman, Omid Bozorg-Haddad, and Xuefeng Chu. "Gradient-based optimizer: A new metaheuristic optimization algorithm." Information Sciences 540 (2020): 131-159.](https://doi.org/10.1016/j.ins.2020.06.037](https://www.sciencedirect.com/science/article/pii/S0020025520306241)))

> Ahmadianfar, I., Bozorg-Haddad, O., & Chu, X. (2020). Gradient-based optimizer: A new metaheuristic optimization algorithm. *Information Sciences*, 540, 131–159. [https://doi.org/10.1016/j.ins.2020.06.037](https://doi.org/10.1016/j.ins.2020.06.037) ([Gradient-Based Optimizer - File Exchange - MATLAB Central](https://www.mathworks.com/matlabcentral/fileexchange/131588-gradient-based-optimizer))
> Ahmadianfar, I., Heidari, A. A., Gandomi, A. H., Chu, X., & Chen, H. (2021). RUN beyond the metaphor: An efficient optimization algorithm based on Runge Kutta method. *Expert Systems with Applications, 181*, 115079. https://doi.org/10.1016/j.eswa.2021.115079([https://www.sciencedirect.com/science/article/pii/S0957417421005200])

## 🔧 Problem Description

We aim to determine optimal **release decisions** for a reservoir system to minimize squared water supply deficits while preventing over-releasing and managing evaporation, storage limits, and possible spill.

- **Input Data:**
  - `Nile` inflow time series (built-in R dataset)
  - Synthetic demands: 80% of inflow + noise
  - Evaporation: 20% of inflow + noise

- **Constraints:**
  - Storage cannot exceed capacity extremes 
  - Releases cannot exceed demand
  - Spill occurs if storage exceeds capacity
---

## ⚙️ Initialization & algorithms Settings

```r
inflow <- as.numeric(Nile)   # Nile R. Annual flow
n_years <- length(inflow)
demand <- rnorm(n_years, inflow * 0.8, inflow * 0.8 * 0.2) # Create synthetic demand
evaporation <- rnorm(n_years, inflow * 0.2, inflow * 0.2 * 0.2) # Create synthetic evaporation

nP <- 100       # Population size
MaxIt <- 200   # Max iterations
lb <- rep(0,n_years)         # Lower bound
max_capacity <- 400 # Maximum storage capacity
min_capacity <- 100
ub <- demand  # Upper bound

```

---

## 🧠 Cost Function

The cost function considers:
- Squared deficits (unmet demand)
- Penalty for over-releasing more than demand

```r
# Returns cost and full system states
evaluate_policy <- function(release)
{
    storage <- numeric(n_years + 1)
    spill <- numeric(n_years)
    deficit <- numeric(n_years)
    storage[1] <- max_capacity / 2
    penalty<-0

    for (t in 1:n_years)
    {
        net_storage = storage[t] + inflow[t] - release[t] - evaporation[t]
        spill[t] <- max(0, net_storage - max_capacity)
        storage[t + 1] <- min(max(net_storage, min_capacity), max_capacity)
        if(net_storage < min_capacity) penalty<-penalty+(min_capacity-net_storage)*100
        deficit[t] <- max(0, demand[t] - release[t])
    }

    total_deficit <- sum(deficit^2)

    return(list(
      cost = total_deficit + sum(penalty),
      storage = storage,
      spill = spill,
      release = release,
      deficit = deficit
    ))
}
```

---

## 🚀 Optimization Execution

```r

obj_function <- function(x) evaluate_policy(x)$cost

source("https://raw.githubusercontent.com/rarabzad/GBO/refs/heads/main/GBO.R")
result_GBO <- GBO(nP, MaxIt, lb, ub, obj_function)
source("https://raw.githubusercontent.com/rarabzad/GBO/refs/heads/main/RUN.R")
result_RUN <- RUN(nP, MaxIt, lb, ub, obj_function)

state_GBO <- evaluate_policy(result_GBO$Best_X)
state_RUN <- evaluate_policy(result_RUN$Best_X)
```

---

## 📈 Visualizations with `ggplot2`

```r
library(ggplot2)
library(patchwork)

# Prepare data
df_main <- data.frame(
  Year = 1:n_years,
  Release = state_GBO$release,
  Demand = demand,
  Storage = state_GBO$storage[-1],
  Spill = state_GBO$spill
)

df_convergence <- data.frame(
  Iteration = 1:length(result_GBO$Convergence_curve),
  BestCost_GBO = result_GBO$Convergence_curve,
  BestCost_RUN = result_RUN$Convergence_curve
)

# 1. Convergence
p1 <- ggplot(df_convergence, aes(x = Iteration)) +
  geom_line(aes(y = BestCost_GBO, color = "GBO"), linewidth = 1) +
  geom_line(aes(y = BestCost_RUN, color = "RUN"), linetype = "dashed", linewidth = 1) +
  labs(title = "Convergence Curve", x = "Iteration", y = "Best Cost") +
  scale_color_manual(name = "Algorithm", values = c("GBO" = "blue", "RUN" = "red")) +
  theme_minimal()

# 2. Release vs Demand
p2 <- ggplot(df_main, aes(x = Year)) +
  geom_line(aes(y = Release), color = "blue") +
  geom_line(aes(y = Demand), color = "red", linetype = "dashed") +
  labs(title = "Release vs Demand", y = "Volume") +
  theme_minimal()

# 3. Storage
p3 <- ggplot(df_main, aes(x = Year, y = Storage)) +
  geom_line(color = "darkgreen") +
  labs(title = "Storage Over Time", y = "Storage (MCM)") +
  theme_minimal()

# 4. Spill
p4 <- ggplot(df_main, aes(x = Year, y = Spill)) +
  geom_col(fill = "darkred") +
  labs(title = "Spill Over Time", y = "Spill (MCM)") +
  theme_minimal()

# Combine all
p1 / p2 / p3 / p4 + plot_layout(ncol = 1)
```

---

## 📌 Notes

- This implementation helps visualize how water is managed over time.
- Optimization balances supply, demand, and operational limits.
- Results show the effectiveness of GBO in water resources planning.

---

## 🧾 Citation

Ahmadianfar, Iman, Omid Bozorg-Haddad, and Xuefeng Chu. "Gradient-based optimizer: A new metaheuristic optimization algorithm." Information Sciences 540 (2020): 131-159.
