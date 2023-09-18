
using Random

# Hull-White parameters
a = 0.1
σ = 0.01

# Swap parameters
fixed_rate = 0.05
maturity = 5.0  # years
dt = 0.01
n_steps = Int(maturity / dt)
n_simulations = 10000

# θ(t) function
θ(t) = 0.03 + 0.02 * t

# Zero coupon bond price under Hull-White
function bond_price(t, T, r, a, σ, θ)
    B = (1 - exp(-a * (T - t))) / a
    A = exp((θ(t) - σ^2 / (2 * a^2)) * (B - (T - t)) - σ^2 / (4 * a) * B^2)
    return A * exp(-B * r)
end

# Monte Carlo simulation for swap valuation
function simulate_swap_value(a, σ, θ, fixed_rate, maturity, dt, n_simulations)
    values = zeros(n_simulations)
    for i in 1:n_simulations
        r = 0.03  # Initial short rate
        value = 0.0
        for t in 0:dt:(maturity-dt)
            # Compute the expected floating payment
            floating_payment = r * dt
            # Compute the fixed payment
            fixed_payment = fixed_rate * dt
            # Compute the net payment
            net_payment = fixed_payment - floating_payment
            # Discount the net payment
            value += net_payment * bond_price(t, maturity, r, a, σ, θ)
            # Update the short rate using Euler discretization
            r += (θ(t) - a * r) * dt + σ * sqrt(dt) * randn()
        end
        values[i] = value
    end
    return mean(values)
end

swap_value = simulate_swap_value(a, σ, θ, fixed_rate, maturity, dt, n_simulations)
println("Swap Value: ", swap_value)
