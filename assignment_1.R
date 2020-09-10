

infrared <- read.table("./infrared.dat", header = TRUE)

f12 <- infrared$F12
logf12 <- log(f12)

d1 <- density(logf12, kernel=c("epanechnikov"))
d2 <- density(logf12, bw = 1/2, kernel = "epanechnikov")
d3 <- density(logf12, bw = .1, kernel = "epanechnikov")

epen_kern <- function(x, h, m = 512) {
	rg = range(x)
	xx <- seq(rg[1] - 0 * h, rg[2] + 0 * h, length.out = m)
	y <- numeric(m)
	
	const = length(x) * h
	y <- sapply(xx, function(z) {sum((abs((z-x)/h) <= 1) * (3/4 * (1 - ((z -x)/(h))^2))) / const}) 
	list(x = xx, y = y)
}

d1
h1 <- epen_kern(logf12, 1)


# Bandwidth selection 
# The Epanechnikov kernel is given by
# K(x) = 3/4(1 - x^2) for |x| in [-1, 1] and 0 else

## Plugin assuming gaussian density and using emperical std var 
K2 <- integrate(function(x) ((3/4) * (1 - x^2))^2,  -1, 1 )
sigma2 <- var(logf12)
hn <- function(n, sigma2, sigmak) {
	((8 * sqrt(pi) * K2$value) / (3 * sigmak^2))^(1/5) * (1/n^(1/5)) * sqrt(sigma2)
}

hn(10, var(logf12))

plot(1:10000, hn(1:10000, var(logf12)))


n <- length(logf12)
h_hat <- hn(n, sd(logf12), 1/5)

# Density comparison
d4 <- density(logf12, bw = h_hat, kernel=c("epanechnikov"))
plot(epen_kern(logf12, h_hat), type="l")
lines(d4, col="red")

## Plugin assuming gaussian density and use silverman var 
sigma_tilde <- min(sd(logf12), IQR(logf12) / 1.34)

hn_2 <- function(n, sigma2) (4/3)^(1/5) * sqrt(sigma2) * n^(-1/5)
h_hat2 <- hn_2(n, sigma_tilde)

# Density comparison
d4 <- density(logf12, bw = h_hat2, kernel=c("epanechnikov"))
plot(epen_kern(logf12, h_hat2), type="l")
lines(d4, col="red")



