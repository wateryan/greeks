/// Black-Scholes formula parameters
/// s0 - underlying price (USD per share)
/// x - strike price (USD per share)
/// sigma - volatility (%)
/// r - continuously compounded risk-free interest rate (%)
/// q - continously compounded dividend yield (%)
/// t - time to expiration (% of year)
use std::f64::consts::E;
use std::f64::consts::PI;
use stats::cnd;

pub fn delta_call(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let cnd = cnd(d1);
    let e = E.powf(-(q * t));
    return e * cnd;
}

pub fn delta_put(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let cnd = cnd(d1);
    let e = E.powf(-(q * t));
    return e * (cnd - 1.0);
}

pub fn gamma(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
    let arg1 = E.powf(-(q * t)) / (s0 * sigma * (t.sqrt()));
    let arg2 = 1.0 / ((2.0 * PI).sqrt());
    let arg3 = E.powf((-d1).powf(2.0)) / 2.0;
    return arg1 * arg2 * arg3;
}

pub fn d1(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let ln = (s0 / x).ln();
    let t_num = t * (r - q + (sigma.powf(2f64) / 2f64));
    return (ln + t_num) / (sigma * t.sqrt());
}

pub fn d2(s0: f64, x: f64, t: f64, r: f64, q: f64, sigma: f64) -> f64 {
    let d1 = d1(s0, x, t, r, q, sigma);
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
    const TIME_TO_EXPIRY: f64 = 23.0 / 365.0;

    const E_D1: f64 = 0.0214;
    const E_D2: f64 = -0.1053;
    const E_CALL_DELTA: f64 = 0.5079;
    const E_PUT_DELTA: f64 = -0.4908;
    const E_GAMMA: f64 = 0.0243;

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
}
