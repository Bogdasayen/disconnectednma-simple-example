
b.data <- list()

# Data on connected RCTs

# Number of studies
b.data$ns <- 5

# Number of treatemnts
b.data$nt <- 5

# Number of events in the 5 connected RCTs (4 have 2 arms, 1 has 3 arms)
b.data$r <- matrix(NA, nrow = 5, ncol = 3)
b.data$r[1,] <- c(94, 92, NA)
b.data$r[2,] <- c(98, 96, NA)
b.data$r[3,] <- c(98, 96, NA)
b.data$r[4,] <- c(95, 97, NA)
b.data$r[5,] <- c(94, 98, 96)
# Number of patients on each arm
b.data$n <- matrix(NA, nrow = 5, ncol =3)
b.data$n[1,] <- c(100, 100, NA)
b.data$n[2,] <- c(100, 100, NA)
b.data$n[3,] <- c(100, 100, NA)
b.data$n[4,] <- c(100, 100, NA)
b.data$n[5,] <- c(100, 100, 100)

# Treatments of each arm of each of the 5 RCTs
b.data$t <- matrix(NA, nrow = 5, ncol =3)
b.data$t[1,] <- c(1, 2, NA)
b.data$t[2,] <- c(1, 2, NA)
b.data$t[3,] <- c(1, 3, NA)
b.data$t[4,] <- c(1, 3, NA)
b.data$t[5,] <- c(1, 2, 3)

# Covariate value on each arm of each RCTRs
b.data$x <- matrix(NA, nrow = 5, ncol =3)
b.data$x[1,] <- c(-0.1336573, -0.2878004, NA)
b.data$x[2,] <- c(0.9109553, -0.1664376, NA)
b.data$x[3,] <- c(0.6923338,  0.2929722, NA)
b.data$x[4,] <- c(1.3306355,  1.1880444, NA)
b.data$x[5,] <- c(1.3141482,  1.2933084, 1.700588)

# Number of arms in 5 RCTs of connected network
b.data$na <- c(2, 2, 2, 2, 3)

# Data for reference/baseline treatments

# Number of arms in connected network with reference treatment
# In this case all 5 RCTs included the reference treatment
b.data$ns.base <- 5

# Number of events in arms on reference treatment
b.data$r.base <- c(94, 98, 98, 95, 94)

# Number of patients in arms on reference treatment
b.data$n.base <- c(100, 100, 100, 100, 100)

# Covariate values in arms on reference treatment
# Corresponds to first column of x[] matrix in connected RCTs
b.data$x.base <- c(-0.1336573, 0.9109553, 0.6923338, 1.3306355, 1.3141482)

# Average value of covariate in arms on reference treatment
b.data$x.base.mean <- c(0.8228831)

# Data for single-arm studies to be included using ALM or reference prediction
# Number of single-arm studies
b.data$ns.single <- 10

# Number of events on single-arm studies
b.data$r.single <- c(100, 96, 98, 87, 99, 91, 36, 78, 99, 83)

# Number of patients on single-arm studies
b.data$n.single <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100)

# Treatments of single-arm studies
b.data$t.single <- c(4, 4, 4, 4, 4, 5, 5, 5, 5, 5)

# Covariate values for single-arm studies
b.data$x.single <- c(-0.43956105, -0.36811619, 0.71098340, -0.13441955, 0.26357280, 2.25849534,
                     0.10642150, -0.06012546, 1.29700372, 0.63297998)

# Data for disconnected RCTs to include using either ALM or reference prediction
# Number of disconnected RCTs
b.data$ns.disc <- 5

# Number of events on arms of disconnected RCTs
b.data$r.disc <- matrix(NA, nrow = 5, ncol = 2)
b.data$r.disc[1,] <- c(99,  100)
b.data$r.disc[2,] <- c(96,  100)
b.data$r.disc[3,] <- c(90,   78)
b.data$r.disc[4,] <- c(98,   90)
b.data$r.disc[5,] <- c(96,   95)

# Number of patients on arms of disconnected RCTs
b.data$n.disc <- matrix(NA, nrow = 5, ncol = 2)
b.data$n.disc[1,] <- c(100, 100)
b.data$n.disc[2,] <- c(100, 100)
b.data$n.disc[3,] <- c(100, 100)
b.data$n.disc[4,] <- c(100, 100)
b.data$n.disc[5,] <- c(100, 100)

# Treatments for arms of disconnected RCTs
b.data$t.disc <- matrix(NA, nrow = 5, ncol = 3)
b.data$t.disc[1,] <- c(4, 5, NA)
b.data$t.disc[2,] <- c(4, 5, NA)
b.data$t.disc[3,] <- c(4, 5, NA)
b.data$t.disc[4,] <- c(4, 5, NA)
b.data$t.disc[5,] <- c(4, 5, NA)

# Covariate values for arms of disconnected RCTs
b.data$x.disc <- matrix(NA, nrow = 5, ncol = 2)
b.data$x.disc[1,] <- c(2.2867433, 1.20999082)
b.data$x.disc[2,] <- c(1.5475021, 1.01273989)
b.data$x.disc[3,] <- c(1.1266772, 0.23397848)
b.data$x.disc[4,] <- c(0.2853316, 0.08022455)
b.data$x.disc[5,] <- c(1.8597435, 0.31436530)

# Number of arms in disconnected RCTs
b.data$na.disc <- c(2, 2, 2, 2, 2)