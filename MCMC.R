model <- BMMmodel(x, k = 2, initialValues = list(S0 = 2), priors = list(kind = "independence", parameter = "priorsFish", hierarchical = "tau"))
control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"), burn.in = 1000, n.iter = 5000, seed = 10)

z <- JAGSrun(x, model = model, control = control, tmp = FALSE, cleanup = TRUE)

plot(z)