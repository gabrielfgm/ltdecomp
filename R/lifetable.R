
lifetable <- function(age, Dx, Px) {
  # Death rates
  mx <- Dx / Px
  # Number of age groups
  ngroup <- length(Dx)

  # Length of each age group
  n <- c(1, 4, rep(5, ngroup - 2))

  # Average person years lived by those dying in the interval
  ax <- n/2
  ax[1] <- 0.35
  ax[ngroup] <- 1/mx[ngroup]

  # STEP 1: From rates to probabilities
  qx <- n*mx / (1 + (n-ax)*mx)
  qx[ngroup] <- 1

  # STEP 2: Survival curve
  # General formula
  lx <- c(1, cumprod(1 - qx[-ngroup]))

  # STEP 3: Life table deaths
  dx <- lx*qx

  # STEP 4: Person years
  Lx <- n * (lx - dx) + ax * dx
  Lx[ngroup] <- lx[ngroup] / mx[ngroup]

  # STEP 5: Person-years above age x
  Tx <- rev(cumsum(rev(Lx)))

  # STEP 6: Life expectancy
  ex <- Tx/lx

  # Create the life table
  LT <- data.frame(age = Age, mx = mx, ax = ax, qx = qx,
                   lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex)

  LT

}




