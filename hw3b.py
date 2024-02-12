import math

def gamma_function(alpha):
    """
    This will calculate the gamma function.
    :param alpha: This is the parameter of the gamma function
    :return: The solution to the gamma function.
    """
    def inte(t):
        return math.exp(-t) * t**(alpha - 1)

    a = 0
    b = 1000  # Can adjust the upper limit
    n = 10000  # Can adjust the number of intervals

    h = (b - a) / n
    result = (inte(a) + inte(b))

    for i in range(1, n, 2):
        result += 1.5 * inte(a + i * h)

    for i in range(2, n - 1, 2):
        result += 2 * inte(a + i * h)

    return (h / 3) * result

def tdp(x, m):
    """
    This calculates the right hand side of the t-distribution probability equation.
    :param x: This is the upper limit of the integral.
    :param m: This is the degrees of freedom.
    :return: The probability value.
    """
    alpha = 1/2 * m + 1/2
    K_m = gamma_function(alpha) / (math.sqrt(m * math.pi) * gamma_function(1/2 * m))

    a = -1000  # Can adjust the lower limit
    n = 10000  # Can adjust the number of intervals

    h = (x - a) / n
    result = (1 + a**2/m)**(-(m + 1)/2)

    for i in range(1, n, 2):
        result += 4 * (1 + (a + i * h)**2/m)**(-(m + 1)/2)

    for i in range(2, n - 1, 2):
        result += 2 * (1 + (a + i * h)**2/m)**(-(m + 1)/2)

    return K_m * (h / 3) * result

def main():
    """
    Once ran, enter the degrees of freedom and 3 desired z values
    this will then calculate and print the probabilities(F(Z)).
    :return: F(z)
    """

    degrees_of_freedom = int(input("Enter the degrees of freedom:"))
    z_values = [float(input(f"Enter z value {i+1} (float): ")) for i in range(3)]

    for z in z_values:
        probability = tdp(z, degrees_of_freedom)
        print(f"For z={z} and {degrees_of_freedom} degrees of freedom, F(z) = {probability}")

if __name__ == "__main__":
    main()