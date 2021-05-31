#----------------------
# Shay Fletcher     318727641
# Nika Tatikishvili 321735433
https://github.com/shayfletcherz/NumericalAnalysisProject4.git
#----------------------

import math
import sympy as sp
from sympy.utilities.lambdify import lambdify

#פונקציה להדפסת הפונקציה המקורית ונגזרתה
def printDerived(function):
    x = sp.symbols('x')
    print("f(x) : ", function)
    derivedFunction = sp.diff(function, x)
    print("f'(x) : ", derivedFunction)
    secondDerived = sp.diff(derivedFunction)
    print("f''(x) : ", secondDerived)

#חישוב נגזרת
def calcDerived(function):
    # calc the derivative from func -> lambdify
    x = sp.symbols('x')
    f_prime = function.diff(x)
    return f_prime

#הפונקציה תחזיר את המספר איטרציות לשיטת החצייה
def calcError(startPoint, endPoint, epsilon):
    return -(math.log(epsilon / (endPoint - startPoint)) / math.log(2))

#שיטת החצייה - תפריד את הטווחים לטווחים קטנים יותר ותחזיר את השורשים
def partition(polynom, startPoint, endPoint, method, choice, epsilon):
    x = sp.symbols('x')
    f = lambdify(x, polynom)
    polynomTag = calcDerived(polynom)
    fTag = lambdify(x, polynomTag)
    i = startPoint
    while i < endPoint: #ריצה על הטווחים
        j = i + 0.1
        if choice is 1:  #ביסקשן
            error = calcError(startPoint, endPoint, epsilon)
            c, iteration = method(f, i, j, error, epsilon)
        elif choice is 2:  #ניוטון-רפסון
            c, iteration = method(polynom, i, j, epsilon)
        else:  #שיטת המיתר
            c, iteration = method(f, i, j, epsilon)

        if c is not None:
            print("root: " + str(c) + " , and " + str(iteration) + " iterations")

        #פעולות על שורש
        if choice is 1:  #ביסקשן
            error = calcError(startPoint, endPoint, epsilon)
            c, iteration = method(fTag, i, j, error, epsilon)  # derived
        elif choice is 2:  #ניוטון-רפסון
            c, iteration = method(polynomTag, i, j, epsilon)
        else:  #שיטת המיתר
            c, iteration = method(fTag, i, j, epsilon)

        if c is not None:
            if -epsilon < f(c) < epsilon:
                print("root: " + str(c) + " , and " + str(iteration) + " iterations")
        i += 0.1

#פונקציית הדפסה בעזרת Bisection
def bisectionMethodPrint(polynom, startPoint, endPoint, choice, epsilon):
    print("Bisection_Method\n")
    partition(polynom, startPoint, endPoint, runBisection, choice, epsilon)

#הפונקציה תציג את כמות השורשים בטווח הנתון
def runBisection(function, startPoint, endPoint, error, epsilon):
    fXl = function(startPoint)
    fXr = function(endPoint)
    if (fXl * fXr) > 0:
        return None, None
    i = -1
    c = startPoint
    while endPoint - startPoint > epsilon and i < error:
        i += 1
        c = (startPoint + endPoint) / 2 #חישוב נקודת אמצע
        if (function(startPoint) * function(c)) > 0:
            startPoint = c
        else:
            endPoint = c
    if i + 1 > error: #במידה ואין שורשים
        print("Could not find root.")
        return None, None
    return c, i

#פונקציית הדפסה בעזרת ניוטון רפסון
def newtonRaphsonMethodPrint(polynom, startPoint, endPoint, choice, epsilon):
    print("Newton Raphson\n")
    partition(polynom, startPoint, endPoint, runNewtonRephson, choice, epsilon)

#הפונקציה תציג את כמות השורשים בטווח הנתון
def runNewtonRephson(f, startA, endB, epsilon):
    x = sp.symbols('x')
    fTag = lambdify(x, calcDerived(f))
    f = lambdify(x, f)
    fXl = f(startA)
    fXr = f(endB)
    if (fXl * fXr) > 0:
        return None, None
    i = -1
    c = (startA + endB) / 2 #חישוב נקודת אמצע
    newC = abs(startA + endB)
    while i < 100:
        i += 1
        temp = newC
        newC = c - (f(c) / fTag(c))
        if abs(newC - c) < epsilon: #במידה ואין שורשים
            return newC, i
        c = temp
    print("Could not find root.")
    return None, None

#פונקציית הדפסה בעזרת שיטת המיתר
def secantMethodPrint(polynom, startPoint, endPoint, choice, epsilon):
    print("Secant method\n")
    partition(polynom, startPoint, endPoint, runSecant, choice, epsilon)

#הפונקציה תציג את כמות השורשים בטווח הנתון
def runSecant(f, startA, endB, epsilon):
    fXl = f(startA)
    fXr = f(endB)
    if (fXl * fXr) > 0:
        return None, None
    c = startA
    newC = endB
    i = -1
    while i < 100:
        i += 1
        temp = newC
        newC = (c * f(newC) - newC * f(c)) / (f(newC) - f(c))
        c = temp
        if abs(newC - c) < epsilon: #במידה ואין שורשים
            return newC, i
    print("Could not find root.")
    return None, None

def Main():
    x = sp.symbols('x')
    f =4 * x ** 3 - 48 * x + 5
    startPoint = 3
    endPoint = 4
    epsilon = 0.0001

    printDerived(f)
    print("\n")
    choice = input("Enter 1 for Bisection method, 2 for Newton Raphson, 3 for Secant method:\n")
    if choice is '1':
        bisectionMethodPrint(f, startPoint, endPoint, 1, epsilon)
    elif choice is '2':
        newtonRaphsonMethodPrint(f, startPoint, endPoint, 2, epsilon)
    elif choice is '3':
        secantMethodPrint(f, startPoint, endPoint, 3, epsilon)
    else:
        bisectionMethodPrint(f, startPoint, endPoint, 1, epsilon)
        newtonRaphsonMethodPrint(f, startPoint, endPoint, 2, epsilon)
        secantMethodPrint(f, startPoint, endPoint, 3, epsilon)

Main()  #הרצה
