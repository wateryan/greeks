/// Black-Scholes formula parameters
/// s0 - underlying price (USD per share)
/// x - strike price (USD per share)
/// sigma - volatility (%)
/// r - continuously compounded risk-free interest rate (%)
/// q - continously compounded dividend yield (%)
/// t - time to expiration (% of year)
/// days_per_year - number of days per year (generally 365)
use std::f64::consts::E;
use std::f64::consts::PI;
use stats::cnd;

/// Calculates the delta of a call option.
///
/// Delta measures the rate of the theoretical option value with respect to the changes in the underlying asset's price.
pub fn delta_call(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let cnd = cnd(d1);
    let e = E.powf(-(q * t));
    return e * cnd;
}

/// Calculates the delta of a put options
///
/// Delta measures the rate of the theoretical option value with respect to the changes in the underlying asset's price.
pub fn delta_put(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let cnd = cnd(d1);
    let e = E.powf(-(q * t));
    return e * (cnd - 1.0);
}

/// Calculates the Gamma for an option
///
/// Gamma measures the rate of change in the delta with respect to the change in the underlying price.
pub fn gamma(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    return gamma_d1(s0, t, q, sigma, d1);
}

pub fn gamma_d1(s0: f64, t: f64, q: f64, sigma: f64, d1: f64) -> f64 {
    let arg1 = E.powf(-(q * t)) / (s0 * sigma * (t.sqrt()));
    let arg2 = one_over_sqrt_pi();
    let arg3 = E.powf((-d1).powf(2.0)) / 2.0;
    return arg1 * arg2 * arg3;
}

/// Calculates the Theta of a call option
///
/// Theta measures the sensitivity of the value of the derivative to the passage of time.
pub fn theta_call(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64, days_per_year: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let arg1 = theta_arg_1(s0, t, q, sigma, d1);
    let d2 = d2_d1(t, sigma, d1);
    let arg2 = theta_arg_2(x, t, r, d2);
    let arg3 = theta_arg_3(s0, t, q, d1);
    return (1.0 / days_per_year) * (arg1 - arg2 + arg3);
}

/// Calculates the Theta of a put option
///
/// Theta measures the sensitivity of the value of the derivative to the passage of time.
pub fn theta_put(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64, days_per_year: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let arg1 = theta_arg_1(s0, t, q, sigma, d1);
    let d2 = d2_d1(t, sigma, d1);
    let arg2 = theta_arg_2(x, t, r, -d2); // d2 is negative for a put
    let arg3 = theta_arg_3(s0, t, q, -d1); // d1 is negative for a put
    return (1.0 / days_per_year) * (arg1 + arg2 - arg3);
}

fn theta_arg_1(s0: f64, t: f64, q: f64, sigma: f64, d1: f64) -> f64 {
    return -(((s0 * sigma * E.powf(-q * t)) / (2.0 * t.sqrt())) * one_over_sqrt_pi() *
             E.powf((-d1.powf(2.0)) / 2.0));
}

fn theta_arg_2(x: f64, t: f64, r: f64, d2: f64) -> f64 {
    return r * x * E.powf(-r * t) * cnd(d2);
}

fn theta_arg_3(s0: f64, t: f64, q: f64, d1: f64) -> f64 {
    return q * s0 * E.powf(-q * t) * cnd(d1);
}

/// Calculates the Vega of a given option
///
/// Vega measures the sensitivity to volatility. Vega is the derivative of the option value with respect to the volatility of the underlying asset.
pub fn vega(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    return vega_d1(s0, t, q, d1);
}

pub fn vega_d1(s0: f64, t: f64, q: f64, d1: f64) -> f64 {
    let mult1 = (1.0 / 100.0) * s0 * E.powf(-(q * t)) * t.sqrt();
    let mult2 = one_over_sqrt_pi();
    let mult3 = E.powf((-d1.powf(2.0) / 2.0));
    return mult1 * mult2 * mult3;
}

/// Calculates the Rho of a call option
///
/// Rho measures the sensitivity to the interest rate. Rho is the derivative of the option value with respect to the risk free interest rate.
pub fn rho_call(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d2_cnd = cnd(d2(s0, x, t, r, q, sigma));
    return (1.0 / 100.0) * x * t * E.powf(-r * t) * d2_cnd;
}

/// Calculates the Rho of a put option
///
/// Rho measures the sensitivity to the interest rate. Rho is the derivative of the option value with respect to the risk free interest rate.
pub fn rho_put(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let neg_d2_cnd = cnd(-d2(s0, x, t, r, q, sigma));
    return -(1.0 / 100.0) * x * t * E.powf(-r * t) * neg_d2_cnd;
}

fn one_over_sqrt_pi() -> f64 {
    return 1.0 / (2.0 * PI).sqrt();
}

pub fn d1(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let ln = (s0 / x).ln();
    let t_num = t * (r - q + (sigma.powf(2f64) / 2f64));
    return (ln + t_num) / (sigma * t.sqrt());
}

// TODO Add overload for providing d1
pub fn d2(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    return d1 - (t.sqrt() * sigma);
}

pub fn d2_d1(t: f64, sigma: f64, d1: f64) -> f64 {
    return d1 - (t.sqrt() * sigma);
}

#[cfg(test)]
mod tests {

    use greeks::*;

    const UNDERLYING: f64 = 64.68;
    const STRIKE: f64 = 65.00;
    const VOL: f64 = 0.5051;
    const INTEREST_RATE: f64 = 0.0150;
    const DIV_YIELD: f64 = 0.0210;
    const DAYS_PER_YEAR: f64 = 365.0;
    const TIME_TO_EXPIRY: f64 = 23.0 / DAYS_PER_YEAR;

    const E_D1: f64 = 0.0214;
    const E_D2: f64 = -0.1053;
    const E_CALL_DELTA: f64 = 0.5079;
    const E_PUT_DELTA: f64 = -0.4908;
    const E_GAMMA: f64 = 0.0243;
    const E_THETA_CALL: f64 = -0.0703;
    const E_THETA_PUT: f64 = -0.0714;
    const E_VEGA: f64 = 0.0647;
    const E_RHO_CALL: f64 = 0.0187;
    const E_RHO_PUT: f64 = -0.0222;

    #[test]
    fn test_d1() {
        let d1 = d1(UNDERLYING,
                    STRIKE,
                    TIME_TO_EXPIRY,
                    INTEREST_RATE,
                    DIV_YIELD,
                    VOL);
        println!("{}", d1);
        let abs = (d1 - E_D1).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_d2() {
        let d2 = d2(UNDERLYING,
                    STRIKE,
                    TIME_TO_EXPIRY,
                    INTEREST_RATE,
                    DIV_YIELD,
                    VOL);
        let abs = (d2 - E_D2).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_delta_call() {
        let call_delta = delta_call(UNDERLYING,
                                    STRIKE,
                                    TIME_TO_EXPIRY,
                                    INTEREST_RATE,
                                    DIV_YIELD,
                                    VOL);
        let abs = (call_delta - E_CALL_DELTA).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_delta_put() {
        let put_delta = delta_put(UNDERLYING,
                                  STRIKE,
                                  TIME_TO_EXPIRY,
                                  INTEREST_RATE,
                                  DIV_YIELD,
                                  VOL);
        let abs = (put_delta - E_PUT_DELTA).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_gamma() {
        let gamma = gamma(UNDERLYING,
                          STRIKE,
                          TIME_TO_EXPIRY,
                          INTEREST_RATE,
                          DIV_YIELD,
                          VOL);
        let abs = (gamma - E_GAMMA).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_theta_call() {
        let theta_call = theta_call(UNDERLYING,
                                    STRIKE,
                                    TIME_TO_EXPIRY,
                                    INTEREST_RATE,
                                    DIV_YIELD,
                                    VOL,
                                    DAYS_PER_YEAR);
        let abs = (theta_call - E_THETA_CALL).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_theta_put() {
        let theta_put = theta_put(UNDERLYING,
                                  STRIKE,
                                  TIME_TO_EXPIRY,
                                  INTEREST_RATE,
                                  DIV_YIELD,
                                  VOL,
                                  DAYS_PER_YEAR);
        println!("{}", theta_put);
        let abs = (theta_put - E_THETA_PUT).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_vega() {
        let vega = vega(UNDERLYING,
                        STRIKE,
                        TIME_TO_EXPIRY,
                        INTEREST_RATE,
                        DIV_YIELD,
                        VOL);
        let abs = (vega - E_VEGA).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_rho_call() {
        let rho_call = rho_call(UNDERLYING,
                                STRIKE,
                                TIME_TO_EXPIRY,
                                INTEREST_RATE,
                                DIV_YIELD,
                                VOL);
        let abs = (rho_call - E_RHO_CALL).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_rho_put() {
        let rho_put = rho_put(UNDERLYING,
                              STRIKE,
                              TIME_TO_EXPIRY,
                              INTEREST_RATE,
                              DIV_YIELD,
                              VOL);
        let abs = (rho_put - E_RHO_PUT).abs();
        assert!(abs < 0.001);
    }

}
