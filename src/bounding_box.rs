use crate::geometry::primitives::FPCoordinate;

#[derive(Debug, Eq, PartialEq)]
pub struct BoundingBox {
    min: FPCoordinate,
    max: FPCoordinate,
}

impl BoundingBox {
    pub fn from_coordinates(coordinates: &[FPCoordinate]) -> Self {
        let mut min_coordinate = FPCoordinate::max();
        let mut max_coordinate = FPCoordinate::min();

        coordinates.iter().for_each(|coordinate| {
            min_coordinate.lat = min_coordinate.lat.min(coordinate.lat);
            min_coordinate.lon = min_coordinate.lon.min(coordinate.lon);
            max_coordinate.lat = max_coordinate.lat.max(coordinate.lat);
            max_coordinate.lon = max_coordinate.lon.max(coordinate.lon);
        });

        BoundingBox {
            min: min_coordinate,
            max: max_coordinate,
        }
    }

    pub fn center(&self) -> FPCoordinate {
        debug_assert!(self.min.lat <= self.max.lat);
        debug_assert!(self.min.lon <= self.max.lon);

        let lat_diff = self.max.lat - self.min.lat;
        let lon_diff = self.max.lon - self.min.lon;

        FPCoordinate {
            lat: self.min.lat + lat_diff / 2,
            lon: self.min.lon + lon_diff / 2,
        }
    }
}

impl From<BoundingBox> for geojson::Bbox {
    fn from(bbox: BoundingBox) -> geojson::Bbox {
        let result = vec![
            bbox.min.lon as f64 / 1000000.,
            bbox.min.lat as f64 / 1000000.,
            bbox.max.lon as f64 / 1000000.,
            bbox.max.lat as f64 / 1000000.,
        ];
        result
    }
}

#[cfg(test)]
pub mod tests {
    use crate::{bounding_box::BoundingBox, geometry::primitives::FPCoordinate};

    #[test]
    pub fn grid() {
        let mut coordinates: Vec<FPCoordinate> = Vec::new();
        for i in 0..100 {
            coordinates.push(FPCoordinate::new(i / 10, i % 10));
        }

        let expected = BoundingBox {
            min: FPCoordinate::new(0, 0),
            max: FPCoordinate::new(9, 9),
        };
        let result = BoundingBox::from_coordinates(&coordinates);
        assert_eq!(expected, result);
    }

    #[test]
    pub fn center() {
        let center = BoundingBox {
            min: FPCoordinate::new_from_lat_lon(33.406637, -115.000801),
            max: FPCoordinate::new_from_lat_lon(33.424732, -114.905286),
        }
        .center();
        assert_eq!(center, FPCoordinate::new(33415684, -114953044));
    }

    #[test]
    pub fn center_with_rounding() {
        let center = BoundingBox {
            min: FPCoordinate::new(0, 0),
            max: FPCoordinate::new(9, 9),
        }
        .center();
        assert_eq!(center, FPCoordinate::new(4, 4));
    }

    #[test]
    pub fn center_without_rounding() {
        let center = BoundingBox {
            min: FPCoordinate::new(0, 0),
            max: FPCoordinate::new(100, 100),
        }
        .center();
        assert_eq!(center, FPCoordinate::new(50, 50));
    }
}
