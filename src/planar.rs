//! Planar stratified media

use crate::material;


#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanarLayer {
    pub height: f64,
    pub material: material::Material
}

impl PlanarLayer {
    pub fn new() -> Self {
        PlanarLayer { height: 0.0, material: material::Material::VACUUM }
    }

    pub fn with_material(height: f64, material: material::Material) -> Self {
        PlanarLayer { height, material }
    }

}

impl Default for PlanarLayer {
    fn default() -> Self {
        Self::new()
    }
}

pub struct PlanarStratifiedMedia<const N: usize> {
    pub layers: [PlanarLayer;N],
}

impl<const N: usize> PlanarStratifiedMedia<N> {
    pub fn new() -> Self {
        PlanarStratifiedMedia { layers: [PlanarLayer::default();N] }
    }
}

impl<const N: usize> Default for PlanarStratifiedMedia<N> {
    fn default() -> Self {
        Self::new()
    }
}
