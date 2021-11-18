from smoothing import smoothing

# Linear regression
print(smoothing.linear_regression_weighted([1,2,3],[3,4,5],[1,1,1]))
# should give slope 1.0, offset 2

data = [0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0]
degree = 4
m = 8
print(smoothing.ms_smoother(data, 0, degree, m))
