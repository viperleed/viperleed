from smoothing import smoothing

# Test linear regression
print("Test linear regression")
print("Result:")
print(smoothing.linear_regression_weighted([1,2,3],[3,4,5],[1,1,1])) # should give slope 1.0, offset 2
print("Expected result:")
print("[2.0, 1.0]]")

print("\n"*2)

# Testing data by Michael
print("Test MS data smoothing")
data = [0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0]
degree = 4
m = 8
print("Input data:")
print(data)
print("Result:")
print(smoothing.ms_smoother(data, 0, degree, m))
print("Expected result:")
print("[-8.31549200e-03  1.38124109e-01  1.83790092e-02 -4.79257597e-02\n"
      " 1.06123277e-01  1.83319786e-01 -3.98482097e-01 -9.58140997e-01\n"
      " 8.05549125e-01  5.18276957e+00  8.36415390e+00  7.23942767e+00\n"
      " 3.60794710e+00  8.63521172e-01 -3.14624208e-01]")
