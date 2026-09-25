# Copyright (C) 2018, 2019, 2020 Dominic O'Kane


# Allow this example to run directly from its category folder.


import time
from financepy.models.process_simulator import VasicekNumericalSchemeTypes
from financepy.models.process_simulator import CIRNumericalSchemeTypes
from financepy.models.process_simulator import HestonNumericalSchemeTypes
from financepy.models.process_simulator import ProcessTypes
from financepy.models.process_simulator import ProcessSimulator

########################################################################################


########################################################################################

# ============================================================================
# 1. FIN PROCESS SIMULATOR
# ============================================================================
# What this section demonstrates:
# Runs the original FinancePy calculation with explicit inputs so the numerical result and the effect of the chosen assumptions can be inspected.

print("\n" + "=" * 78)
print("1. FIN PROCESS SIMULATOR")
print("=" * 78)


# ============================================================================
# FINANCEPY EXAMPLES - ProcessSimulator
# ============================================================================

num_paths = 20000
num_annual_steps = 100
seed = 1919
t = 1.0
model_sim = ProcessSimulator()
print_paths = False

print("######################## GBM NORMAL ###############################")
sigma = 0.10
stock_price = 100.0
drift = 0.04
scheme = GBMNumericalSchemeTypes.NORMAL
model_params = (stock_price, drift, sigma, scheme)
start = time.time()
paths = model_sim.get_process(ProcessTypes.GBM_PROCESS, t, model_params, num_annual_steps, num_paths, seed)
end = time.time()
elapsed = end - start
print("PROCESS", "TIME")
print("GBM NORMAL", elapsed)
if print_paths:
    print(paths)

print("######################## GBM ANTITHETIC ###########################")
sigma = 0.10
stock_price = 100.0
drift = 0.04
scheme = GBMNumericalSchemeTypes.ANTITHETIC
model_params = (stock_price, drift, sigma, scheme)
start = time.time()
paths = model_sim.get_process(ProcessTypes.GBM_PROCESS, t, model_params, num_annual_steps, num_paths, seed)
end = time.time()
elapsed = end - start
print("GBM ANTITHETIC", elapsed)
if print_paths:
    print(paths)

print("###################### HESTON euler ###############################")
stock_price = 100.0
v0 = 0.05
kappa = 0.50
theta = 0.05
sigma = 0.90
rho = -0.9
scheme = HestonNumericalSchemeTypes.EULER
model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
start = time.time()
paths = model_sim.get_process(ProcessTypes.HESTON_PROCESS, t, model_params, num_annual_steps, num_paths, seed)
end = time.time()
elapsed = end - start
print("HESTON EULER", elapsed)
if print_paths:
    print(paths)

print("###################### HESTON EULERLOG ############################")
stock_price = 100.0
v0 = 0.05
kappa = 0.50
theta = 0.05
sigma = 0.90
rho = -0.9
scheme = HestonNumericalSchemeTypes.EULERLOG
model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
start = time.time()
paths = model_sim.get_process(ProcessTypes.HESTON_PROCESS, t, model_params, num_annual_steps, num_paths, seed)
end = time.time()
elapsed = end - start
print("HESTON EULERLOG", elapsed)
if print_paths:
    print(paths)

print("###################### HESTON QUADEXP #############################")
stock_price = 100.0
v0 = 0.05
kappa = 0.50
theta = 0.05
sigma = 0.90
rho = -0.9
scheme = HestonNumericalSchemeTypes.QUADEXP
model_params = (stock_price, drift, v0, kappa, theta, sigma, rho, scheme)
start = time.time()
paths = model_sim.get_process(ProcessTypes.HESTON_PROCESS, t, model_params, num_annual_steps, num_paths, seed)
end = time.time()
elapsed = end - start
print("HESTON QUADEXP", elapsed)
if print_paths:
    print(paths)

print("######################## VASICEK NORMAL ###########################")
r0 = 0.05
kappa = 0.50
theta = 0.05
sigma = 0.90
scheme = VasicekNumericalSchemeTypes.NORMAL
model_params = (r0, kappa, theta, sigma, scheme)
start = time.time()
paths = model_sim.get_process(
    ProcessTypes.VASICEK_PROCESS,
    t,
    model_params,
    num_annual_steps,
    num_paths,
    seed,
)
end = time.time()
elapsed = end - start
print("VASICEK_NORMAL", elapsed)
if print_paths:
    print(paths)

print("####################### VASICEK ANTITHETIC ########################")
r0 = 0.05
kappa = 0.50
theta = 0.05
sigma = 0.90
scheme = VasicekNumericalSchemeTypes.ANTITHETIC
model_params = (r0, kappa, theta, sigma, scheme)
start = time.time()
paths = model_sim.get_process(
    ProcessTypes.VASICEK_PROCESS,
    t,
    model_params,
    num_annual_steps,
    num_paths,
    seed,
)
end = time.time()
elapsed = end - start
print("VASICEK_NORMAL ANTI", elapsed)
if print_paths:
    print(paths)

print("############################# CIR #################################")
r0 = 0.05
kappa = 0.50
theta = 0.05
sigma = 0.90
scheme = CIRNumericalSchemeTypes.MILSTEIN
model_params = (r0, kappa, theta, sigma, scheme)
start = time.time()
paths = model_sim.get_process(ProcessTypes.CIR_PROCESS, t, model_params, num_annual_steps, num_paths, seed)
end = time.time()
elapsed = end - start
print("CIR", elapsed)
if print_paths:
    print(paths)
