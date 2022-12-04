# biquad-filter

Biquad IIR filters for real time audio.
Based on the [filter design by Robert Bristow-Johnson](https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt).

## Example
```Rust
use biquad_filter::Biquad;

let mut filter = Biquad::default();
let mut filter = Biquad::<f32>::default();
filter.prepare(44100).unwrap();
filter.set(FilterType::Peak, 100., 2., 1.).unwrap();

let input = vec![0.; 100];
let mut output = vec![0.; 100];
filter.process(&input, &mut output);
```