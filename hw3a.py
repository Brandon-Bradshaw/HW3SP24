import math
def SymPosDef(matrix):
    """
    This checks to see if the Matrix is symmetrical and positive definite.
    :param matrix: This is the Matrix that is going to be checked.
    :return: True if the Matrix is symmetrical and positive definite, False if not.
    """
    rows, cols = len(matrix), len(matrix[0])
    for i in range(rows):
        for j in range(cols):
            if matrix[i][j] != matrix[j][i]:
                return False
    return True
    row, cols = len(matrix), len(matrix[0])
    for i in range(rows):
        minor = [row[:i + 1] for row in matrix[:i + 1]]
        determinant = 1
        for j in range(i + 1):
            determinant *= minor[j][j]
        if determinant <= 0:
            return False
    return True
def ChoDec(matrix):
    """
    This uses Cholesky Decomposition on a symmetric positive definite matrix.
    :param matrix: The matrix that is being checked.
    :return: L and x. L is the lower triangular matrix and x is the solution vector.
    """
    n = len(matrix)
    L = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            if i == j:
                sum_value = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = math.sqrt(matrix[i][i] - sum_value)
            else:
                sum_value = sum(L[i][k] * L[j][k] for k in range (j))
                L[i][j] = (matrix[i][j] - sum_value) / L[j][j]
    return L
def SolveEqCho(matrix, vector):
    """
    This uses the Cholesky method to solve the system of linear equations.
    :param matrix: The matrix that is being used.
    :param vector: The right hand side.
    :return: The solution vector.
    """

    n = len(matrix)
    L = ChoDec(matrix)
    y = [0] * n
    for i in range(n):
        y[i] = vector[i] - sum(L[i][j] * y[j] for j in range(i))
    x = [0] * n
    for i in reversed(range(n)):
        x[i] = (y[i] - sum(L[j][i] * x[j] for j in range(i + 1, n))) / L[i][i]
    return x, "Cholesky Method"
def SolveEq(matrix, vector):
    """
    This uses the Doolittle method to solve the system of linear equations.
    :param matrix: The matrix that is being checked.
    :param vector: The right hand side.
    :return: The solution vector.
    """
    n = len(matrix)
    for k in range (n - 1):
        for i in range(k + 1, n):
            factor = matrix[i][k] / matrix [k][k]
            matrix[i][k] = factor
            for j in range (k +1, n):
                matrix[i][k] -= factor * matrix[k][j]
    y = [0] * n
    for i in range(n):
        y[i] = vector[i] - sum(matrix[i][j] * y[j] for j in range(i))
    x = [0] * n
    for i in reversed(range(n)):
        x[i] = (y[i] - sum(matrix[i][j] * x[j] for j in range (i + 1, n))) / matrix[i][i]
    return x, "Doolittle Method"
matrix1 = [
    [1, -1, 3, 2],
    [-1, 5, -5, -2],
    [3, -5, 19, 3],
    [2, -2, 3,21]
]
vector1 = [15, -35, 94, 1]
matrix2 = [
    [4, 2, 4, 0],
    [2, 2, 3, 2],
    [4, 3, 6, 3],
    [0, 2, 3, 9]
]
vector2 = [20, 36, 60, 122]

if SolveEqCho(matrix1, vector1):
    solution, method_used = SolveEqCho(matrix1, vector1)
    solution2, method_used2 = SolveEqCho(matrix2, vector2)
else:
    solution, method_used = SolveEq(matrix1, vector1)
    solution2, method_used2 = SolveEq(matrix2, vector2)

print("Solution vector1:", solution)
print("Solution vector2:", solution2)
print("Method used1:", method_used)
print("Method used2:", method_used2)


