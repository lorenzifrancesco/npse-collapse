import numpy as np
from scipy.optimize import newton
from functools import partial
import matplotlib.pyplot as plt


def find_non_monotone_index(arr):
    if len(arr) < 2:
        return -1  # Array is too short to determine monotonicity

    increasing = None

    for i in range(1, len(arr)):
        if arr[i] == arr[i - 1]:
            continue  # Skip equal elements

        if increasing is None:
            increasing = arr[i] > arr[i - 1]
        elif (increasing and arr[i] < arr[i - 1]) or (not increasing and arr[i] > arr[i - 1]):
            return i

    return -1  # Array is monotone


GAMMAC = 0.67

# Define the function and its derivative


def f(mu, gamma):
    return (1-mu)**(3/2) - 3/2 * (1-mu)**(1/2) + 3/(2*np.sqrt(2)) * gamma


def wf(spi, z, mu):
    # spi = np.abs(spi)
    print(spi)
    return 1/np.sqrt(2) * 1/np.sqrt(1-mu) * np.arctan(np.sqrt((spi - mu)/(1-mu))) - 1/np.sqrt(2) * 1/np.sqrt(1+mu) * np.arctanh(np.sqrt((spi - mu)/(1+mu))) - z


# Find the root using SciPy's newton method
gammas = np.linspace(0, 0.8, 100)
mus = np.zeros_like(gammas)
first_failure = 0
for i in range(len(gammas)):
    fus = partial(f, gamma=gammas[i])
    try:
        mus[i] = newton(fus, 1.0)
    except:
        if first_failure == 0:
            first_failure = i
            print("============", first_failure)
        mus[i] = 0.0
    # print(mus[i])

ic = find_non_monotone_index(mus)
print(f"Gamma collapse: {gammas[ic]:4.3e}")
print(f"Mu    collapse: {mus[ic]:4.3e}")
plt.figure(figsize=(5, 4))
plt.plot(gammas, mus)
plt.plot(gammas, 1-gammas**2/2)
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$\mu$")
plt.grid()
plt.tight_layout()
plt.savefig("media/mus.pdf")

space = np.linspace(-5, 5, 1000)
sigma = np.zeros_like(space)

for i in range(len(space)):
    wfz = partial(wf, z=space[i], mu=0.8)
    try:
        sigma[i] = newton(wfz, 0.1)
    except:
        if first_failure == 0:
            first_failure = i
            print("============", first_failure)
        sigma[i] = 0.0
    # print(mus[i])

plt.figure(figsize=(5, 4))
plt.plot(space, sigma)
plt.xlabel(r"$\zeta$")
plt.ylabel(r"$\sigma$")
plt.grid()
plt.tight_layout()
plt.savefig("media/sigma.pdf")

# middlec = partial(middle, mu = mus[ic])
# spic = newton(middlec, 100)
