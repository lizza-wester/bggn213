student1 <- c(100, 100, 100, 100, 100, 100, 100, 90)
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
student3 <- c(90, NA, NA, NA, NA, NA, NA, NA)

grade1 <- function(x) {
  low = which.min(x)
  score_minus_low = sum(x) - low
  score1 = score_minus_low/length(x)
  print(score1)
}