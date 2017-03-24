# greeks

[![Build Status](https://travis-ci.org/wateryan/greeks.svg?branch=master)](https://travis-ci.org/wateryan/greeks)
[![Crates.io Status](https://img.shields.io/crates/v/greeks.svg)](https://crates.io/crates/greeks)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This project is a simple implementation of options pricing and greeks for Rust.

Goals:

* Be fast
* Be accurrate
* Not have any external dependencies.

## Supported Calculations

### Greeks

#### First Order

* Delta
* Lambda
* Rho
* Theta
* Vega

#### Second Order

* Gamma

### Pricing

* European call option
* European put option

### Valution

* Call option at expiry
* Put option at expiry