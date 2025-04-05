use integrate::prelude::legendre_rule;
use nannou::glam::Vec2;

pub trait Trajectory: Send + Sync {
    fn position(&self, t: f32, theta: f32) -> Vec2;
    fn velocity(&self, t: f32, theta: f32) -> Vec2 {
        let h = 1e-4;
        let dx = (self.position(t + h, theta).x - self.position(t, theta).x) / h;
        let dy = (self.position(t + h, theta).y - self.position(t, theta).y) / h;
        Vec2::new(dx, dy)
    }
    fn a(t: f32, theta: f32) -> f32;

    fn s(&self, t: f32, theta: f32) -> f32 {
        let ds = |u: f32| {
            let vel = self.velocity(u, theta);
            (vel.x.powi(2) + vel.y.powi(2)).powf(0.5)
        };
        let upper_limit: f32 = t;
        static LOWER_LIMIT: f32 = 0.;
        static N: usize = 50;
        legendre_rule(ds, LOWER_LIMIT, upper_limit, N) as f32
    }

    fn n(&self, t: f32, theta: f32) -> Vec2 {
        let vel = self.velocity(t, theta);
        let ds = (vel.x.powi(2) + vel.y.powi(2)).powf(0.5);
        if ds < 1e-6 {
            Vec2::new(1., 0.)
        } else {
            Vec2::new(-vel.y, vel.x) / ds
        }
    }
}

pub trait WaveFunction: Trajectory {
    fn normal_offset(&self, t: f32, theta: f32, ta: f32) -> f32 {
        let current_s = self.s(t, theta);
        Self::a(t, theta) * (current_s / ta).sin()
    }

    fn offset(&self, t: f32, theta: f32, ta: f32) -> Vec2 {
        let normal = self.n(t, theta);
        self.position(t, theta) + normal * self.normal_offset(t, theta, ta)
    }
}

impl<T: Trajectory> WaveFunction for T {}
