#' Construct a lifetable
#'
#' @param age A vector of ages or age ranges, e.g. '0-1', '1-4', '5-9', etc
#' @param Dx A numeric vector of death counts for the age groups.
#' @param Px A numeric vector of population sizes for the age groups.
#'
#' @return A dataframe containing a completed life-table.
#' The column 'age' returns the age groups.
#' Column 'mx' is the mortality rate at age x.
#' Column 'ax' is the estimated mortality rate for individuals that die within the age range of x, e.g. individuals that die between '1-4' for age '1-4'.
#' Column 'qx' is the probability of dying at age x.
#' Column 'lx' is the magnitude of the population at x. It has a fixed radix of 1.
#' Column 'dx' is the number of lifetable deaths.
#' Column 'Lx' is the number of person years lived in age bin x.
#' Column 'Tx' is the number of person years lived above age x.
#' Column 'ex' is life-expectancy at age x.
#' @export
#'
#' @examples
#' # Age groups
#' Age <-  c("0", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34",
#' "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69",
#' "70-74", "75-79", "80-84", "85-89", "90-94", "95+")
#'
#' # Death counts
#' Dx <- c(110934, 17693, 9168, 7378, 12192, 13356, 14223, 19210, 29136, 42961, 64310,
#'         90630, 116759, 153504, 196682, 223796, 220065, 185305, 120414, 50299, 13887)
#'
#' # mid year population
#' Px <- c(4126560, 16195304, 18659141, 16815965, 13287434, 10803165, 10870165,
#'         11951709, 12508316, 11567216, 10528878, 9696502, 8595947, 7111897, 6186763,
#'         4661136, 2977347, 1518206, 648581, 170653, 44551)
#'
#' lifetable(age = Age, Dx = Dx, Px = Px)
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
  LT <- data.frame(age = age, mx = mx, ax = ax, qx = qx,
                   lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex)

  LT

}




