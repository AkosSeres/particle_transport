use eframe::egui::{self, RichText};

use crate::photon::set_detector;

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
    pub simulation_running: bool,
}

impl MyApp {
    pub fn init(&self) {
        set_detector(
            self.arguments.radius,
            self.arguments.height,
            self.arguments.density,
        );
    }

    fn start_stop_simulation(&mut self) {
        if !self.simulation_running {
            set_detector(
                self.arguments.radius,
                self.arguments.height,
                self.arguments.density,
            );
        }

        self.simulation_running = !self.simulation_running;
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
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
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.radius, 0.5..=12.0).text("cm"),
                    );
                    ui.end_row();
                    ui.label("Height: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.height, 0.5..=12.0).text("cm"),
                    );
                    ui.end_row();
                });

                ui.heading(RichText::from("Emitter position coordinates:").size(18.0));
                egui::Grid::new("emitter_position_grid").show(ui, |ui| {
                    ui.label("X: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.rx, -10.0..=10.0).text("cm"),
                    );
                    ui.end_row();
                    ui.label("Y: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.ry, -10.0..=10.0).text("cm"),
                    );
                    ui.end_row();
                    ui.label("Z: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.rz, -10.0..=10.0).text("cm"),
                    );
                    ui.end_row();
                });

                ui.heading(RichText::from("Other settings:").size(18.0));
                egui::Grid::new("other_settings_gui").show(ui, |ui| {
                    ui.label("Photon energy: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.energy, 200.0..=2500.0).text("keV"),
                    );
                    ui.end_row();
                    ui.label("FWHM: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.fwhm, 0.0..=100.0).text("keV"),
                    );
                    ui.end_row();
                    ui.label("Detector density: ");
                    ui.add_enabled(
                        !self.simulation_running,
                        egui::Slider::new(&mut self.arguments.density, 0.1..=20.0).text("g/cm³"),
                    );
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

                ui.vertical_centered(|ui| {
                    ui.add_space(15.0);

                    let start_stop_text = if self.simulation_running {
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
                    if self.simulation_running {
                        ui.add_space(5.0);
                        ui.add(egui::Spinner::default());
                    }
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let sin = (0..1000).map(|i| {
                let x = i as f64 * 0.01;
                egui::plot::Value::new(x, x.sin())
            });
            let line = egui::plot::Line::new(egui::plot::Values::from_values_iter(sin));
            egui::plot::Plot::new("main_plot_spectrum")
                .allow_zoom(true)
                .show(ui, |plot_ui| plot_ui.line(line));
        });
    }
}

impl Default for MyApp {
    fn default() -> Self {
        Self {
            arguments: MyArgs::default(),
            simulation_running: false,
        }
    }
}
