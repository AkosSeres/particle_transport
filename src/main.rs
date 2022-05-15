#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use crate::photon::set_detector;
use clap::Parser;
use eframe::egui::{self, RichText};

mod photon;
pub mod rand_gen;
mod vec3;

/// Simulates the transport of the photons of a monoenergetic gamma-source inside a scintillation detector.
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Radius of the detector in cm
    #[clap(long, default_value_t = 2.5)]
    radius: f64,

    /// Height of the detector in cm
    #[clap(short, long, default_value_t = 3.0)]
    height: f64,

    /// Energy of the emitted photons in keV
    #[clap(short, long, default_value_t = 661.7)]
    energy: f64,

    /// The density of the detector material, in g/cm³
    #[clap(short, long, default_value_t = 3.67)]
    density: f64,

    /// Full width at half maximum (FWHM) of the detector, in keV
    #[clap(long, default_value_t = 6.0)]
    fwhm: f64,

    /// X coordinate of the emitter
    #[clap(long, default_value_t = 3.0)]
    rx: f64,
    /// Y coordinate of the emitter
    #[clap(long, default_value_t = -3.0)]
    ry: f64,
    /// Z coordinate of the emitter
    #[clap(long, default_value_t = 2.0)]
    rz: f64,

    /// Number of particles to evaluate
    #[clap(short, long, default_value_t = 10_000)]
    particlenum: u64,
}

fn main() {
    let args = Args::parse();
    set_detector(args.radius, args.height, args.density);

    let mut options = eframe::NativeOptions::default();
    options.vsync = true;
    options.transparent = false;
    options.maximized = true;
    eframe::run_native(
        "Photon transport",
        options,
        Box::new(|_cc| Box::new(MyApp { arguments: args })),
    );
}

struct MyApp {
    arguments: Args,
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
                    ui.add(egui::Slider::new(&mut self.arguments.radius, 0.5..=12.0).text("cm"));
                    ui.end_row();
                    ui.label("Height: ");
                    ui.add(egui::Slider::new(&mut self.arguments.height, 0.5..=12.0).text("cm"));
                    ui.end_row();
                });

                ui.heading(RichText::from("Emitter position coordinates:").size(18.0));
                egui::Grid::new("emitter_position_grid").show(ui, |ui| {
                    ui.label("X: ");
                    ui.add(egui::Slider::new(&mut self.arguments.rx, -10.0..=10.0).text("cm"));
                    ui.end_row();
                    ui.label("Y: ");
                    ui.add(egui::Slider::new(&mut self.arguments.ry, -10.0..=10.0).text("cm"));
                    ui.end_row();
                    ui.label("Z: ");
                    ui.add(egui::Slider::new(&mut self.arguments.rz, -10.0..=10.0).text("cm"));
                    ui.end_row();
                });

                ui.heading(RichText::from("Other settings:").size(18.0));
                egui::Grid::new("other_settings_gui").show(ui, |ui| {
                    ui.label("Photon energy: ");
                    ui.add(
                        egui::Slider::new(&mut self.arguments.energy, 200.0..=2500.0).text("keV"),
                    );
                    ui.end_row();
                    ui.label("FWHM: ");
                    ui.add(egui::Slider::new(&mut self.arguments.fwhm, 0.0..=100.0).text("keV"));
                    ui.end_row();
                    ui.label("Detector density: ");
                    ui.add(
                        egui::Slider::new(&mut self.arguments.density, 0.1..=20.0).text("g/cm³"),
                    );
                    ui.end_row();
                });

                egui::ScrollArea::vertical().show(ui, |ui| {
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
                });

                ui.vertical_centered(|ui| {
                    ui.add_space(15.0);

                    let start_button = egui::widgets::Button::new(
                        RichText::new("Start simulation").strong().size(18.0),
                    );
                    if ui.add(start_button).clicked() {
                        println!("Simulation is wanting to start!")
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
