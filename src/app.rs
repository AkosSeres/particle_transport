use crate::{
    photon::{set_default_energy, set_detector, Photon, F},
    photon_emitter::PhotonEmitter,
    rand_gen::RandGen,
    vec3::Vector,
};
use eframe::{
    egui::{self, RichText},
    epaint::Color32,
};
use std::sync::{
    atomic::{
        AtomicBool, AtomicU64,
        Ordering::{Relaxed, SeqCst},
    },
    Arc,
};
#[cfg(not(target_arch = "wasm32"))]
use std::thread;

/// Simulates the transport of the photons of a monoenergetic gamma-source inside a scintillation detector.
#[derive(Debug)]
pub struct MyArgs {
    /// Radius of the detector in cm
    pub radius: f64,
    /// Height of the detector in cm
    pub height: f64,
    /// Energy of the emitted photons in keV
    pub energy: f64,
    /// The density of the detector material, in g/cm³
    pub density: f64,
    /// Full width at half maximum (FWHM) of the detector, in keV
    pub fwhm: f64,
    /// X coordinate of the emitter
    pub rx: f64,
    /// Y coordinate of the emitter
    pub ry: f64,
    /// Z coordinate of the emitter
    pub rz: f64,
}

impl Default for MyArgs {
    fn default() -> Self {
        Self {
            radius: 2.5,
            height: 3.0,
            energy: 661.7,
            density: 3.67,
            fwhm: 6.0,
            rx: 3.0,
            ry: -3.0,
            rz: 2.0,
        }
    }
}

pub struct MyApp {
    pub arguments: MyArgs,
    pub simulation_running: Arc<AtomicBool>,
    /// The number of channels of the spectrometer
    channel_number: usize,
    channels: Arc<Vec<AtomicU64>>,
    start_instant: Option<instant::Instant>,
    end_instant: Option<instant::Instant>,
    logscale: bool,
}

impl MyApp {
    pub fn init(&self) {
        set_detector(
            self.arguments.radius,
            self.arguments.height,
            self.arguments.density,
        );
        set_default_energy(self.arguments.energy);
    }

    fn get_max_energy(&self) -> f64 {
        self.arguments.energy + self.arguments.fwhm * 10.0
    }

    fn start_stop_simulation(&mut self) {
        if !self.simulation_running.load(SeqCst) {
            set_detector(
                self.arguments.radius,
                self.arguments.height,
                self.arguments.density,
            );
            set_default_energy(self.arguments.energy);
            self.start_simulation();
        } else {
            self.stop_simulation();
        }
    }

    fn start_simulation(&mut self) {
        self.reset_spectrum();
        let max_energy = self.get_max_energy();
        self.start_instant = Some(instant::Instant::now());
        self.end_instant = None;

        #[cfg(not(target_arch = "wasm32"))]
        for _ in 0..num_cpus::get() {
            let channels = self.channels.clone();
            let simulation_running = self.simulation_running.clone();
            let energy = self.arguments.energy;
            let (rx, ry, rz) = (self.arguments.rx, self.arguments.ry, self.arguments.rz);
            let fwhm = self.arguments.fwhm;
            let photon_emit = PhotonEmitter::from_params(
                self.arguments.radius,
                self.arguments.height,
                rx,
                ry,
                rz,
            );
            let pos_vector = Vector::<F>::new(rx, ry, rz);
            thread::spawn(move || loop {
                for _ in 0..100000 {
                    let mut random_photon = Photon {
                        energy,
                        pos: pos_vector,
                        dir: photon_emit.gen_photon_dir(),
                    };
                    let mut energy_hit_size = random_photon.simulate();
                    if energy_hit_size <= 0.0 {
                        continue;
                    }
                    for _ in 0..12 {
                        energy_hit_size += (F::rand() - 0.5) * fwhm;
                    }
                    let channel_width = max_energy / (channels.len() as f64);
                    let idx = (energy_hit_size / channel_width).floor() as usize;
                    if idx >= channels.len() {
                        continue;
                    }
                    channels[idx].fetch_add(1, Relaxed);
                }

                if !simulation_running.load(Relaxed) {
                    return;
                }
            });
        }

        self.simulation_running.store(true, SeqCst);
    }

    fn stop_simulation(&mut self) {
        self.simulation_running.store(false, SeqCst);
        self.end_instant = Some(instant::Instant::now());
    }

    fn reset_spectrum(&mut self) {
        self.channels.iter().for_each(|ch| ch.store(0, SeqCst));
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        let simulation_running = self.simulation_running.load(Relaxed);
        ctx.set_visuals(egui::style::Visuals::dark());
        egui::SidePanel::right("main_settings_panel")
            .resizable(false)
            .show(ctx, |ui| {
                ui.heading(
                    RichText::from("Particle transport simulation")
                        .size(24.0)
                        .strong(),
                );

                ui.heading(RichText::from("Detector size:").size(18.0));
                egui::Grid::new("detector_size_settings_grid").show(ui, |ui| {
                    ui.label("Radius: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.radius, 0.5..=12.0).text("cm"),
                    );
                    ui.end_row();
                    ui.label("Height: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.height, 0.5..=12.0).text("cm"),
                    );
                    ui.end_row();
                });

                ui.heading(RichText::from("Emitter position coordinates:").size(18.0));
                egui::Grid::new("emitter_position_grid").show(ui, |ui| {
                    ui.label("X: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.rx, -10.0..=10.0).text("cm"),
                    );
                    ui.end_row();
                    ui.label("Y: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.ry, -10.0..=10.0).text("cm"),
                    );
                    ui.end_row();
                    ui.label("Z: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.rz, -10.0..=10.0).text("cm"),
                    );
                    ui.end_row();
                });

                ui.heading(RichText::from("Other settings:").size(18.0));
                egui::Grid::new("other_settings_gui").show(ui, |ui| {
                    ui.label("Photon energy: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.energy, 1.0..=20000.0).text("keV"),
                    );
                    ui.end_row();
                    ui.label("FWHM: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.fwhm, 0.0..=100.0).text("keV"),
                    );
                    ui.end_row();
                    ui.label("Detector density: ");
                    ui.add_enabled(
                        !simulation_running,
                        egui::Slider::new(&mut self.arguments.density, 0.1..=20.0).text("g/cm³"),
                    );
                    ui.end_row();
                    ui.label("Logarithmic scale: ");
                    ui.add(egui::Checkbox::new(&mut self.logscale, ""));
                    ui.end_row();
                });

                /*egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading(RichText::from("Views:").size(18.0));
                    let halfrange = self
                        .arguments
                        .radius
                        .max(self.arguments.height / 2.0)
                        .max(self.arguments.rx.abs())
                        .max(self.arguments.ry.abs())
                        .max(self.arguments.rz.abs());
                    egui::CollapsingHeader::new("Top view").show(ui, |ui| {
                        let detector_circle = (0..51).map(|i| {
                            let x = i as f64 * 2.0 * std::f64::consts::PI / 50.0;
                            egui::plot::Value::new(
                                self.arguments.radius * x.cos(),
                                self.arguments.radius * x.sin(),
                            )
                        });
                        let line = egui::plot::Line::new(egui::plot::Values::from_values_iter(
                            detector_circle,
                        ));
                        let point = egui::plot::Points::new(egui::plot::Values::from_values_iter(
                            [egui::plot::Value::new(self.arguments.rx, self.arguments.ry)]
                                .into_iter(),
                        ))
                        .radius(5.0);
                        egui::plot::Plot::new("top_view_plot")
                            .allow_zoom(true)
                            .view_aspect(1.0)
                            .data_aspect(1.0)
                            .allow_boxed_zoom(false)
                            .allow_drag(false)
                            .allow_scroll(false)
                            .allow_zoom(false)
                            .center_x_axis(true)
                            .center_y_axis(true)
                            .include_x(halfrange)
                            .include_y(halfrange)
                            .show_x(false)
                            .show_y(false)
                            .show(ui, |plot_ui| {
                                plot_ui.line(line);
                                plot_ui.points(point);
                            });
                    });

                    egui::CollapsingHeader::new("Front view").show(ui, |ui| {
                        let r = self.arguments.radius;
                        let h = self.arguments.height;
                        let detector_shape = [
                            egui::plot::Value::new(-r / 2.0, -h / 2.0),
                            egui::plot::Value::new(-r / 2.0, h / 2.0),
                            egui::plot::Value::new(r / 2.0, h / 2.0),
                            egui::plot::Value::new(r / 2.0, -h / 2.0),
                            egui::plot::Value::new(-r / 2.0, -h / 2.0),
                        ]
                        .into_iter();
                        let line = egui::plot::Line::new(egui::plot::Values::from_values_iter(
                            detector_shape,
                        ));
                        let point = egui::plot::Points::new(egui::plot::Values::from_values_iter(
                            [egui::plot::Value::new(self.arguments.rx, self.arguments.rz)]
                                .into_iter(),
                        ))
                        .radius(5.0);
                        egui::plot::Plot::new("front_view_plot")
                            .allow_zoom(true)
                            .view_aspect(1.0)
                            .data_aspect(1.0)
                            .allow_boxed_zoom(false)
                            .allow_drag(false)
                            .allow_scroll(false)
                            .allow_zoom(false)
                            .center_x_axis(true)
                            .center_y_axis(true)
                            .include_x(halfrange)
                            .include_y(halfrange)
                            .show_x(false)
                            .show_y(false)
                            .show(ui, |plot_ui| {
                                plot_ui.line(line);
                                plot_ui.points(point);
                            });
                    });

                    egui::CollapsingHeader::new("Side view").show(ui, |ui| {
                        let r = self.arguments.radius;
                        let h = self.arguments.height;
                        let detector_shape = [
                            egui::plot::Value::new(-r / 2.0, -h / 2.0),
                            egui::plot::Value::new(-r / 2.0, h / 2.0),
                            egui::plot::Value::new(r / 2.0, h / 2.0),
                            egui::plot::Value::new(r / 2.0, -h / 2.0),
                            egui::plot::Value::new(-r / 2.0, -h / 2.0),
                        ]
                        .into_iter();
                        let line = egui::plot::Line::new(egui::plot::Values::from_values_iter(
                            detector_shape,
                        ));
                        let point = egui::plot::Points::new(egui::plot::Values::from_values_iter(
                            [egui::plot::Value::new(self.arguments.ry, self.arguments.rz)]
                                .into_iter(),
                        ))
                        .radius(5.0);
                        egui::plot::Plot::new("front_view_plot")
                            .allow_zoom(true)
                            .view_aspect(1.0)
                            .data_aspect(1.0)
                            .allow_boxed_zoom(false)
                            .allow_drag(false)
                            .allow_scroll(false)
                            .allow_zoom(false)
                            .center_x_axis(true)
                            .center_y_axis(true)
                            .include_x(halfrange)
                            .include_y(halfrange)
                            .show_x(false)
                            .show_y(false)
                            .show(ui, |plot_ui| {
                                plot_ui.line(line);
                                plot_ui.points(point);
                            });
                    });
                });*/

                let hits = self.channels.iter().map(|ch| ch.load(Relaxed)).sum::<u64>();

                ui.vertical_centered(|ui| {
                    ui.add_space(15.0);

                    let start_stop_text = if simulation_running {
                        "Stop simulation"
                    } else {
                        "Start simulation"
                    };
                    let start_button = egui::widgets::Button::new(
                        RichText::new(start_stop_text).strong().size(18.0),
                    );
                    if ui.add(start_button).clicked() {
                        self.start_stop_simulation();
                    }
                    if simulation_running {
                        ui.add_space(5.0);
                        ui.add(egui::Spinner::default());
                    }
                    ui.label(format!("{}", hits));
                });

                if self.start_instant.is_some() {
                    let elpased_in_sec = if self.end_instant.is_some() {
                        self.end_instant
                            .unwrap()
                            .duration_since(self.start_instant.unwrap())
                            .as_secs_f64()
                    } else {
                        instant::Instant::now()
                            .duration_since(self.start_instant.unwrap())
                            .as_secs_f64()
                    };
                    let rate = if elpased_in_sec == 0.0 {
                        0.0
                    } else {
                        hits as f64 / elpased_in_sec
                    };
                    let rate_text = if rate > 1000000.0 {
                        format!("{:3.2}M/s", (rate / 1000000.0))
                    } else if rate > 1000.0 {
                        format!("{:3.2}k/s", (rate / 1000.0))
                    } else {
                        format!("{:3.2}/s", rate.round())
                    };
                    egui::Grid::new("simulation_info").show(ui, |ui| {
                        ui.add(egui::widgets::Label::new(
                            RichText::new(format!("Hitrate: {}", rate_text))
                                .strong()
                                .size(18.0),
                        ));
                    });
                }
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let channel_width = self.get_max_energy() / self.channel_number as f64;

            let bars = if self.logscale {
                self.channels
                    .iter()
                    .enumerate()
                    .map(|(ch_num, ch_count)| {
                        egui::plot::Bar::new(ch_num as f64 * channel_width, {
                            let value = ch_count.load(Relaxed) as f64;
                            if value.is_normal() {
                                value.log10()
                            } else {
                                0.0
                            }
                        })
                        .fill(Color32::BLUE)
                    })
                    .collect()
            } else {
                self.channels
                    .iter()
                    .enumerate()
                    .map(|(ch_num, ch_count)| {
                        egui::plot::Bar::new(
                            ch_num as f64 * channel_width,
                            ch_count.load(Relaxed) as f64,
                        )
                        .fill(Color32::BLUE)
                    })
                    .collect()
            };

            let spectrum_chart = egui::plot::BarChart::new(bars);

            egui::plot::Plot::new("main_plot_spectrum")
                .allow_zoom(true)
                .include_x(self.get_max_energy())
                .show(ui, |plot_ui| plot_ui.bar_chart(spectrum_chart));
        });

        #[cfg(target_arch = "wasm32")]
        if simulation_running {
            let photon_emit = PhotonEmitter::from_params(
                self.arguments.radius,
                self.arguments.height,
                self.arguments.rx,
                self.arguments.ry,
                self.arguments.rz,
            );
            let start_time = instant::Instant::now();
            let max_energy = self.get_max_energy();
            loop {
                for _ in 0..100000 {
                    let mut random_photon = Photon {
                        energy: self.arguments.energy,
                        pos: Vector::<F>::new(
                            self.arguments.rx,
                            self.arguments.ry,
                            self.arguments.rz,
                        ),
                        dir: photon_emit.gen_photon_dir(),
                    };
                    let mut energy_hit_size = random_photon.simulate();
                    if energy_hit_size <= 0.0 {
                        continue;
                    }
                    for _ in 0..12 {
                        energy_hit_size += (F::rand() - 0.5) * self.arguments.fwhm;
                    }
                    let channel_width = max_energy / (self.channels.len() as f64);
                    let idx = (energy_hit_size / channel_width).floor() as usize;
                    if idx >= self.channels.len() {
                        continue;
                    }
                    self.channels[idx].fetch_add(1, Relaxed);
                }
                let end_time = instant::Instant::now();
                if (end_time - start_time).as_micros() > 30_000 {
                    break;
                }
            }
        }
    }
}

impl Default for MyApp {
    fn default() -> Self {
        let default_channel_number = 1024;
        Self {
            arguments: MyArgs::default(),
            simulation_running: Arc::new(AtomicBool::new(false)),
            channel_number: default_channel_number,
            channels: Arc::new(
                (0..default_channel_number)
                    .map(|_| AtomicU64::new(0))
                    .collect(),
            ),
            start_instant: None,
            end_instant: None,
            logscale: false,
        }
    }
}
