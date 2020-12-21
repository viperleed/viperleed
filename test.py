import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import scipy

def optimizerHelper(array, func):
    """Reshapes a 1D array to 2D and passes it on to func"""
    return func(array.reshape(1,-1))

dim = 3
# fake data
xvals = []
y = []
for i in range(0,100):
    x = np.random.rand(dim)*30-15
    xvals.append(x)
    y.append((10 * x[0] ** 2 - 20 * x[0] + 
              3 * x[1] ** 2 + 12 * x[1] + 
              0.03 * x[2] ** 2 +
              np.random.rand()*200-100))
X = np.array(xvals).reshape(-1, dim)
# y = [(5 * x ** 2 - 3 * x + 10 + np.random.rand()*100) for x in X]

polyreg = make_pipeline(
        PolynomialFeatures(degree=2),
        LinearRegression()
        )
polyreg.fit(X, y)
X_min = scipy.optimize.minimize(optimizerHelper, np.zeros(3), 
                                args=(polyreg.predict,),method='L-BFGS-B',
                                bounds=[(-10,10) for i in range(0,dim)])
print(X_min.x)
print(polyreg.score())

# prediction and plot along the two axes
for i in range(dim):
    vals = 500
    X_test = np.tile(X_min.x, (vals,1))
    X_test[:,i] = np.linspace(-20,20,num=vals)
    y_pred = polyreg.predict(X_test)
    curv = (polyreg.named_steps['linearregression'].coef_[
            polyreg.named_steps['polynomialfeatures']
            .get_feature_names().index('x{}^2'.format(i))])

    plt.scatter(X[:,i], y, label="data")
    plt.plot(X_test[:,i], y_pred, label="prediction, min = ({:.2f}, {:.2f})"
             .format(X_test[:,i][np.where(
                        y_pred == np.amin(y_pred))[0][0]], np.amin(y_pred)))

    plt.title("Dimension {}: curvature = {:.3f}".format(i, curv))
    plt.legend(loc='upper left')
    plt.show()

