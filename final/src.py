from math import pi, sin, sqrt
import numpy as np

def gause_seidel_iteration(A: np.array, b: np.array, epsilon: float):
    iter_times = 0
    n = len(A)
    x_1 = np.zeros(n)
    x_2 = np.ones(n)
    while True:
        iter_times += 1
        x_1 = x_2.__copy__()
        for i in range(n):
            s = 0
            for j in range(n):
                s += A[i][j] * x_2[j]
            x_2[i] = (b[i] - s + A[i][i] * x_2[i]) / A[i][i]
        if np.max(np.abs(x_2 - x_1)) < epsilon:
            break
    return x_2, iter_times


def Jacobi_iteration(A:np.array, b:np.array, epsilon: float):
    iter_times = 0
    x_1 = np.zeros(len(A))
    x_2 = np.ones(len(A))
    while True:
        iter_times += 1
        x_1 = x_2.__copy__()
        for i in range(len(A)):
            s = 0
            for j in range(len(A)):
                s += A[i][j] * x_1[j]
            x_2[i] = (b[i] - s + A[i][i] * x_1[i]) / A[i][i]
        if np.max(np.abs(x_2 - x_1)) < epsilon:
            break
    return x_2, iter_times

def run_question_2(n: int, print_result: bool = False):
    def f(x): 
        return pi * pi * sin(pi*x)
    def u_e(x):
        return sin(pi*x)
    A = 2 * np.eye(n-1)
    for i in range(n-2):
        A[i][i+1] = -1
        A[i+1][i] = -1
    b = [f(i/n)/n/n for i in range(1, n)]
    u_h_Jacobi, iter_times_Jacobi = Jacobi_iteration(A, b, 1e-10)
    e_h_Jacobi = sqrt(sum([(u_h_Jacobi[i] - u_e((i+1)/n))**2 for i in range(n-1)]))
    u_h_Gauss, iter_times_Gauss = gause_seidel_iteration(A, b, 1e-10)
    e_h_Gauss = sqrt(sum([(u_h_Gauss[i] - u_e((i+1)/n))**2 for i in range(n-1)]))
    if print_result:
        print(f"n={n}, Jacobi: {u_h_Jacobi}")
        print(f"n={n}, Gauss Seidel: {u_h_Gauss}, Error: {e_h_Gauss}")
    print(f"n={n}, Jacobi Error: {e_h_Jacobi}")
    print(f"n={n}, Gauss Seidel Error: {e_h_Gauss}")
    return (e_h_Jacobi, e_h_Gauss), (iter_times_Jacobi, iter_times_Gauss)


question_2_n = [10, 20, 40, 80, 160]
question_2_errors = []
question_2_iters = []
for n in question_2_n:
    error, iters = run_question_2(n, True)
    question_2_errors.append(error)
    question_2_iters.append(iters)
print("Question 2 error table:", question_2_errors)
print("Question 2 iter times table:", question_2_iters)


def newton_iteration(n:int , epsilon: float, Jacobi_function, F_function):
    X = np.zeros(n-1)
    iter_times = 0
    while True:
        iter_times += 1
        Jacobi_matrix = Jacobi_function(X)
        delta = np.linalg.inv(Jacobi_matrix).dot(F_function(X).reshape(n-1, 1))
        X = X - delta.reshape(n-1)
        if np.max(np.abs(delta)) < epsilon:
            break
    return X, iter_times

def run_question_6(n: int, print_result: bool = False):
    h = 1/n
    def f(x): 
        return pi * pi * sin(pi*x) + sin(pi*x) ** 3
    
    def u_e(x):
        return sin(pi*x)

    def Jacobi(input:np.array):
        res = np.zeros((n-1, n-1))
        for i in range(n-1):
            res[i][i] = 2 + 3 * input[i] * input[i] * h * h
            if i > 0:
                res[i][i-1] = -1
            if i < n-2:
                res[i][i+1] = -1
        return res

    def F(X:np.array):
        res = np.zeros(n-1)
        res[0] = 2 * X[0] - X[1] + X[0] ** 3 * h * h - f(h) * (h ** 2)
        res[n-2] = 2 * X[n-2] - X[n-3] + X[n-2] ** 3 * h * h - f(1 - h) * (h ** 2)
        for i in range(1, n-2):
            res[i] = 2 * X[i] - X[i-1] - X[i+1] + X[i] ** 3 * h * h - f((i+1) * h) * (h ** 2)
        return res
    
    X, iter_times = newton_iteration(n, 1e-8, Jacobi, F)
    if print_result:
        print(f"n={n}, Newton: {X}")
    print(f"n={n}, Newton iter times: {iter_times}")
    error = sqrt(sum([(X[i] - u_e((i+1)/n))**2 for i in range(n-1)]))
    return error, iter_times

question_6_n = [10, 20, 40, 80, 160]
question_6_errors = []
question_6_iters = []
for n in question_6_n:
    error, iters = run_question_6(n, True)
    question_6_errors.append(error)
    question_6_iters.append(iters)
print("Question 6 error table:", question_6_errors)
print("Question 6 iter times table:", question_6_iters)