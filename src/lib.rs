use num_traits::Float;
use std::f64::consts::PI;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum BiquadError {
    #[error("the sample rate must be set first")]
    NoSampleRate,
    #[error("the frequency is higher than nyqist")]
    FrequencyOverNyqist,
    #[error("the frequency is 0 or lower than 0")]
    FrequencyTooLow,
    #[error("q is lower than zero")]
    NegativeQ,
    #[error("fatal number conversion error")]
    Fatal,
}

enum FilterType {
    Lowpass,
    Highpass,
    Bandpass1,
    Bandpass2,
    Notch,
    Allpass,
    Peak,
    Lowshelf,
    Highshelf,
}

#[derive(Default, Debug, Clone)]
struct Coefficients<F: Float> {
    sample_rate: F,
    a0: F,
    a1: F,
    a2: F,
    b0: F,
    b1: F,
    b2: F,
}

impl<F: Float> Coefficients<F> {
    pub fn set_sample_rate(&mut self, sample_rate: u32) -> Result<(), BiquadError> {
        self.sample_rate = F::from(sample_rate).ok_or(BiquadError::Fatal)?;
        Ok(())
    }

    pub fn set(
        &mut self,
        filter_type: FilterType,
        frequency: f64,
        gain_db: f64,
        q: f64,
    ) -> Result<(), BiquadError> {
        if self.sample_rate == F::zero() {
            return Err(BiquadError::NoSampleRate);
        }
        if 2.0 * frequency > self.sample_rate.to_f64().ok_or(BiquadError::Fatal)? {
            return Err(BiquadError::FrequencyOverNyqist);
        }
        if frequency < 1.0 {
            return Err(BiquadError::FrequencyTooLow);
        }
        if q < 0.0 {
            return Err(BiquadError::NegativeQ);
        }

        let a = f64::powf(10., gain_db / 40.);
        let omega =
            2. * PI * frequency as f64 / self.sample_rate.to_f64().ok_or(BiquadError::Fatal)?;
        let sin = f64::sin(omega);
        let cos = f64::cos(omega);
        let alpha = sin / 2. * q;
        let beta = 2.0 * f64::sqrt(a) * alpha;

        match filter_type {
            FilterType::Lowpass => {
                self.b0 = F::from((1. - cos) / 2.).ok_or(BiquadError::Fatal)?;
                self.b1 = F::from(1. - cos).ok_or(BiquadError::Fatal)?;
                self.b2 = F::from((1. - cos) / 2.).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Highpass => {
                self.b0 = F::from((1. + cos) / 2.).ok_or(BiquadError::Fatal)?;
                self.b1 = F::from(-(1. + cos)).ok_or(BiquadError::Fatal)?;
                self.b2 = F::from((1. + cos) / 2.).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Bandpass1 => {
                self.b0 = F::from(q * alpha).ok_or(BiquadError::Fatal)?;
                self.b1 = F::zero();
                self.b2 = F::from(-q * alpha).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Bandpass2 => {
                self.b0 = F::from(alpha).ok_or(BiquadError::Fatal)?;
                self.b1 = F::zero();
                self.b2 = F::from(-alpha).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Notch => {
                self.b0 = F::one();
                self.b1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.b2 = F::one();
                self.a0 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Allpass => {
                self.b0 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
                self.b1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.b2 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from(1. + alpha).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Peak => {
                self.b0 = F::from(1. + alpha * a).ok_or(BiquadError::Fatal)?;
                self.b1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.b2 = F::from(1. - alpha * a).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from(1. + alpha / a).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(-2. * cos).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from(1. - alpha / a).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Lowshelf => {
                self.b0 =
                    F::from(a * ((a + 1.0) - (a - 1.0) * cos + beta)).ok_or(BiquadError::Fatal)?;
                self.b1 =
                    F::from(2.0 * a * ((a - 1.0) - (a + 1.0) * cos)).ok_or(BiquadError::Fatal)?;
                self.b2 =
                    F::from(a * ((a + 1.0) - (a - 1.0) * cos - beta)).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from((a + 1.0) + (a - 1.0) * cos + beta).ok_or(BiquadError::Fatal)?;
                self.a1 =
                    F::from(-2.0 * ((a - 1.0) + (a + 1.0) * cos)).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from((a + 1.0) + (a - 1.0) * cos - beta).ok_or(BiquadError::Fatal)?;
            }
            FilterType::Highshelf => {
                self.b0 =
                    F::from(a * ((a + 1.0) + (a - 1.0) * cos + beta)).ok_or(BiquadError::Fatal)?;
                self.b1 =
                    F::from(-2.0 * a * ((a - 1.0) + (a + 1.0) * cos)).ok_or(BiquadError::Fatal)?;
                self.b2 =
                    F::from(a * ((a + 1.0) + (a - 1.0) * cos - beta)).ok_or(BiquadError::Fatal)?;
                self.a0 = F::from((a + 1.0) - (a - 1.0) * cos + beta).ok_or(BiquadError::Fatal)?;
                self.a1 = F::from(2.0 * ((a - 1.0) - (a + 1.0) * cos)).ok_or(BiquadError::Fatal)?;
                self.a2 = F::from((a + 1.0) - (a - 1.0) * cos - beta).ok_or(BiquadError::Fatal)?;
            }
        }
        Ok(())
    }
}

#[derive(Default)]
struct Biquad<F: Float> {
    coefficients: Coefficients<F>,
    x1: F,
    x2: F,
    y1: F,
    y2: F,
}

impl<F: Float> Biquad<F> {
    pub fn set(
        &mut self,
        filter_type: FilterType,
        frequency: f64,
        gain_db: f64,
        q: f64,
    ) -> Result<(), BiquadError> {
        self.coefficients.set(filter_type, frequency, gain_db, q)
    }

    pub fn prepare(&mut self, sample_rate: u32) -> Result<(), BiquadError> {
        self.coefficients.set_sample_rate(sample_rate)
    }

    pub fn process(&mut self, input: &[F], output: &mut [F]) {
        output
            .iter_mut()
            .zip(input)
            .for_each(|(out_sample, in_sample)| *out_sample = self.tick(*in_sample));
    }

    pub fn reset(&mut self) {
        self.x1 = F::zero();
        self.x2 = F::zero();
        self.y1 = F::zero();
        self.y2 = F::zero();
    }

    #[inline]
    pub fn tick(&mut self, input: F) -> F {
        let out = self.coefficients.b0 * input
            + self.coefficients.b1 * self.x1
            + self.coefficients.b2 * self.x2
            - self.coefficients.a1 * self.y1
            - self.coefficients.a2 * self.y2;

        self.x2 = self.x1;
        self.x1 = input;
        self.y2 = self.y1;
        self.y1 = out;

        out
    }

    pub fn set_coefficients(&mut self, coefficients: Coefficients<F>) {
        self.coefficients = coefficients;
    }

    pub fn coefficients(&self) -> Coefficients<F> {
        self.coefficients.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let mut filter = Biquad::<f32>::default();
        filter.prepare(44100).unwrap();
        filter.set(FilterType::Peak, 100., 2., 1.).unwrap();
    }
}
