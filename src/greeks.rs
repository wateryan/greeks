/// Black-Scholes formula parameters
/// s0 - underlying price (USD per share)
/// x - strike price (USD per share)
/// sigma - volatility (%)
/// r - continuously compounded risk-free interest rate (%)
/// q - continously compounded dividend yield (%)
/// t - time to expiration (% of year)
use std::f64::consts::E;
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

    const UNDERLYING: f64 = 36.07;
    const STRIKE: f64 = 35.00;
    const VOL: f64 = 0.4825;
    const INTEREST_RATE: f64 = 0.01;
    const DIV_YIELD: f64 = 0.00;
    const TIME_TO_EXPIRY: f64 = 26.0 / 365.0;

    #[test]
    fn test_d1() {
        let expected_d1: f64 = 0.3038;
        let d1 = d1(UNDERLYING,
                    STRIKE,
                    TIME_TO_EXPIRY,
                    INTEREST_RATE,
                    DIV_YIELD,
                    VOL);
        let abs = (d1 - expected_d1).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_d2() {
        let expected_d2 = 0.1751;
        let d2 = d2(UNDERLYING,
                    STRIKE,
                    TIME_TO_EXPIRY,
                    INTEREST_RATE,
                    DIV_YIELD,
                    VOL);
        let abs = (d2 - expected_d2).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_delta_call() {
        let expected_delta = 0.6194;
        let call_delta = delta_call(UNDERLYING,
                                    STRIKE,
                                    TIME_TO_EXPIRY,
                                    INTEREST_RATE,
                                    DIV_YIELD,
                                    VOL);
        let abs = (call_delta - expected_delta).abs();
        assert!(abs < 0.001);
    }

    #[test]
    fn test_delta_put() {
        let expected_delta = -0.3806;
        let put_delta = delta_put(UNDERLYING,
                                  STRIKE,
                                  TIME_TO_EXPIRY,
                                  INTEREST_RATE,
                                  DIV_YIELD,
                                  VOL);
        let abs = (put_delta - expected_delta).abs();
        println!("{}", put_delta);
        assert!(abs < 0.001);
    }
}
