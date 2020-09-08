
infrared <- read.table("./infrared.dat", header = TRUE)

# 2.2
# The Epanechnikov kernel is given by
# K(x) = 3/4(1 - x^2) for |x| in [-1, 1] and 0 else


f12 <- infrared$F12
logf12 <- log(f12)

d1 <- density(logf12, kernel=c("epanechnikov"))
d2 <- density(logf12, bw = 1/2, kernel = "epanechnikov")
d3 <- density(logf12, bw = .1, kernel = "epanechnikov")

plot(d1, col="blue")
lines(d2, col="red")
lines(d3, col="green")

# 2.3

epen_kern <- function(x, h, m = 512) {
	rg = range(x)
	xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
	y <- numeric(m)
	
	const = length(x) * h
	y <- sapply(xx, function(z) {sum((abs((z-x)/h) <= 1) * (3/4 * (1 - ((z -x)/(h))^2))) / const
		}) 
	list(x = xx, y = y)
}

h2 <- epen_kern(logf12, h = 1/2)
plot(h2, col="blue", type="l")
lines(d2, col="red")


system.time(density(logf12, bw = 1/2, kernel = "epanechnikov"))
system.time(epen_kern(logf12, h = 1/2))



plot(abs(density(logf12, bw = 1/2, kernel = "epanechnikov")$y - epen_kern(logf12, h = 1/2)$y))
plot(abs(density(logf12, bw = .2, kernel = "epanechnikov")$y - epen_kern(logf12, h = .2)$y))
plot(abs(density(logf12, bw = .1, kernel = "epanechnikov")$y - epen_kern(logf12, h = .1)$y))


# 2.4
library(microbenchmark)
x <- rnorm(2^13)
k <- seq(2^5, 2^13, by=(abs(2^3-2^13))/10)

bench <- microbenchmark(
	density(x[1:k[1]], .2),
	density(x[1:k[2]], .2),
	density(x[1:k[3]], .2),
	density(x[1:k[4]], .2),
	density(x[1:k[5]], .2),
	density(x[1:k[6]], .2),
	density(x[1:k[7]], .2),
	density(x[1:k[8]], .2),
	density(x[1:k[9]], .2),	
	density(x[1:k[10]], .2)
)

bench



