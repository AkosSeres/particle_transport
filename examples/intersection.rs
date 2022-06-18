use gnuplot::AxesCommon;
use particle_transport::photon::{set_default_energy, set_detector, Photon, PhotonLocation};
use particle_transport::vec3::Vector;

fn main() {
    set_detector(1.0, 1.0, 5.0);
    set_default_energy(700.0);
    let photon_emitter = Vector::<f64>::new(-2.0, -2.0, -2.0);
    let mut ps = Vec::<Vector<f64>>::new();
    for _ in 0..100000 {
        let photon = Photon {
            energy: 700.0,
            pos: photon_emitter,
            dir: Vector::<f64>::random_isotropic_normed(),
        };
        let location = photon.intersect_detector();
        match location {
            PhotonLocation::OutsideMisses => (),
            PhotonLocation::OutsideInto(dist) => {
                let intersect = photon.pos + photon.dir * dist;
                ps.push(intersect);
            }
            PhotonLocation::Inside(_) => (),
        }
    }

    // generate (x,y,z) matrices
    let x: Vec<f64> = ps.iter().map(|p| p.x).collect();
    let y: Vec<f64> = ps.iter().map(|p| p.y).collect();
    let z: Vec<f64> = ps.iter().map(|p| p.z).collect();

    use gnuplot::{Caption, Color, Figure};

    //let x = [0u32, 1, 2];
    //let y = [3u32, 4, 5];
    let mut fg = Figure::new();
    fg.axes3d()
        .points(x, y, z, &[Caption(""), Color("black")])
        .set_z_range(gnuplot::Fix(-5.0), gnuplot::Fix(5.0))
        .set_y_range(gnuplot::Fix(-5.0), gnuplot::Fix(5.0))
        .set_x_range(gnuplot::Fix(-5.0), gnuplot::Fix(5.0));
    fg.show().expect("Failed to display the plot.");
}
