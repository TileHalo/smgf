
pub mod material;
pub mod planar;
pub mod telegrapher;



pub type SpectralFn = Box<dyn Fn(num::complex::Complex64) -> num::complex::Complex64>;
