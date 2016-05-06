context("partition maps")

test_that("partition map has correct dimensions", {
    order = 2
    N = 10^order
    set.seed(123)
    x=rnorm(N)
    y=rnorm(N)
    pm = partition.map(x,y)
    expect_true(is.list(pm))
    expect_true(length(pm[,3])==N)
    expect_true(nchar(pm[,3][1])==order+1)  # we add 1 to account for 'p' character
})
