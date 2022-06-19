use crate::{
    photon::{set_default_energy, set_detector, Photon},
    photon_emitter::PhotonEmitter,
    rand_gen::RandGen,
    vec3::Vector,
};
use atomic_float::AtomicF64;
use eframe::{
    egui::{self, RichText},
    epaint::Color32,
};
#[cfg(not(target_arch = "wasm32"))]
use std::thread;
use std::{
    io::Write,
    sync::{
        atomic::{
            AtomicBool, AtomicU64,
            Ordering::{Relaxed, SeqCst},
        },
        Arc,
    },
};

/// Contains the properties of the detector and the emitter.
/// The emitter emits monoenergetic photons (they all have the same energy)
/// and the detector is a cylinder centered at the origin, parallel to the z-axis.
#[derive(Debug)]
pub struct SimulationArgs {
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

impl Default for SimulationArgs {
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

/// Contains data regarding the statistics of the
/// efficiency of the detector.
///
/// The values have to be atomic because they will be updated
/// from multiple threads.
#[derive(Default)]
struct SimulationStatistics {
    intrinsic_efficiency_sum: AtomicF64,
    intrinsic_efficiency_square_sum: AtomicF64,
    total_efficiency_sum: AtomicF64,
    total_efficiency_square_sum: AtomicF64,
    total_photon_count: AtomicU64,
    photon_count_with_detector: AtomicU64,
}

/// This function executes the simulation of individual photons
/// and collects statistical data about the efficiency of the detector.
fn run_simulation_cycles(
    energy: f64,
    pos_vector: Vector<f64>,
    photon_emit: &PhotonEmitter,
    deviation: f64,
    max_energy: f64,
    channels: &Arc<Vec<AtomicU64>>,
    simulation_statistics: &Arc<SimulationStatistics>,
) {
    let cycle_count = 100_000;
    let mut photon_count_with_detector = 0;
    let mut total_efficiency_sum = 0.0;
    let mut total_efficiency_square_sum = 0.0;
    let mut intrinsic_efficiency_sum = 0.0;
    let mut intrinsic_efficiency_square_sum = 0.0;

    // Simulate [cycle_count] photons
    for _ in 0..cycle_count {
        // Generate a random photons using the random emitter
        let mut random_photon = Photon {
            energy,
            pos: pos_vector,
            dir: photon_emit.gen_photon_dir(),
        };
        // Run the simulation
        let sim_results = random_photon.simulate();

        // Calculate efficiencies
        let current_total_efficiency = sim_results.energy_transfered / energy;
        total_efficiency_sum += current_total_efficiency;
        total_efficiency_square_sum += current_total_efficiency * current_total_efficiency;
        if sim_results.hit_detector {
            photon_count_with_detector += 1;
            let current_intrinsic_efficiency = sim_results.energy_transfered / energy;
            intrinsic_efficiency_sum += current_intrinsic_efficiency;
            intrinsic_efficiency_square_sum +=
                current_intrinsic_efficiency * current_intrinsic_efficiency;
        }

        let mut energy_hit_size = sim_results.energy_transfered;
        if energy_hit_size <= 0.0 || !sim_results.hit_detector {
            // Continue if the photon did not expend any energy or if it did not hit the detector
            continue;
        }

        // Append a random amount of energy to the energy hit size,
        // to simulate the detector's resolution
        energy_hit_size += f64::rand_normal_sum12() * deviation;
        if energy_hit_size <= 0.0 {
            continue;
        }

        // Register the energy hit in the corresponding channel
        let channel_width = max_energy / (channels.len() as f64);
        let idx = (energy_hit_size / channel_width).floor() as usize;
        if idx >= channels.len() {
            continue;
        }
        channels[idx].fetch_add(1, Relaxed);
    }

    // Store the recorded statistics into the atomic variables
    simulation_statistics
        .total_photon_count
        .fetch_add(cycle_count, Relaxed);
    simulation_statistics
        .photon_count_with_detector
        .fetch_add(photon_count_with_detector, Relaxed);
    let solid_angle = photon_emit.get_solid_angle();
    simulation_statistics.total_efficiency_sum.fetch_add(
        total_efficiency_sum * (solid_angle / (4.0 * std::f64::consts::PI)),
        Relaxed,
    );
    simulation_statistics.total_efficiency_square_sum.fetch_add(
        total_efficiency_square_sum
            * (solid_angle * solid_angle / (16.0 * std::f64::consts::PI * std::f64::consts::PI)),
        Relaxed,
    );
    simulation_statistics
        .intrinsic_efficiency_sum
        .fetch_add(intrinsic_efficiency_sum, Relaxed);
    simulation_statistics
        .intrinsic_efficiency_square_sum
        .fetch_add(intrinsic_efficiency_square_sum, Relaxed);
}

/// Organizes the application, renders the UI and manages execution of the simulations.
pub struct MyApp {
    /// Simulation arguments
    arguments: SimulationArgs,
    // True, if the simulation is currently running, else false
    simulation_running: Arc<AtomicBool>,
    /// The number of channels of the spectrometer
    channel_number: usize,
    /// The atomic integers storing the number of hits in individual channels
    channels: Arc<Vec<AtomicU64>>,
    /// The instant when the possible simulation started
    start_instant: Option<instant::Instant>,
    /// The instant when the last simulation ended
    end_instant: Option<instant::Instant>,
    /// Whether we want to display the spectrum with a logarithmic scale
    logscale: bool,
    /// The statistical data of the simulation
    simulation_statistics: Arc<SimulationStatistics>,
}

impl MyApp {
    /// Sets the initial parameters of the detector
    pub fn init(&self) {
        set_detector(
            self.arguments.radius,
            self.arguments.height,
            self.arguments.density,
        );
        set_default_energy(self.arguments.energy);
    }

    /// Return the maximum energy of spectrometer channels.
    /// Depends on the maximum energy of the photons and the resolution of the detector (FWHM).
    fn get_max_energy(&self) -> f64 {
        self.arguments.energy + self.arguments.fwhm * 3.0
    }

    /// Converts the FWHM into standard deviation.
    fn get_standard_deviation_from_fwhm(&self) -> f64 {
        self.arguments.fwhm / 2.35482004503
    }

    /// Either starts or stops the simulation depending on the current state.
    fn start_stop_simulation(&mut self) {
        if !self.simulation_running.load(SeqCst) {
            // Reset if starting a new simulation
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

    /// Starts the simulation.
    ///
    /// On native platforms, this function spawns new threads to calculate.
    /// On WASM, this function does not start new threads, but sets some variables
    /// to indicate that the simulation must be run in the main loop, on the main thread.
    fn start_simulation(&mut self) {
        self.reset_spectrum();
        let max_energy = self.get_max_energy();
        self.start_instant = Some(instant::Instant::now());
        self.end_instant = None;

        #[cfg(not(target_arch = "wasm32"))]
        // Start as many threads as the number of logical cores
        for _ in 0..num_cpus::get() {
            let channels = self.channels.clone();
            let simulation_running = self.simulation_running.clone();
            let energy = self.arguments.energy;
            let (rx, ry, rz) = (self.arguments.rx, self.arguments.ry, self.arguments.rz);
            let deviation = self.get_standard_deviation_from_fwhm();
            let photon_emit = PhotonEmitter::from_params(
                self.arguments.radius,
                self.arguments.height,
                rx,
                ry,
                rz,
            );
            let pos_vector = Vector::<f64>::new(rx, ry, rz);
            let simulation_statistics = self.simulation_statistics.clone();
            thread::spawn(move || loop {
                run_simulation_cycles(
                    energy,
                    pos_vector,
                    &photon_emit,
                    deviation,
                    max_energy,
                    &channels,
                    &simulation_statistics,
                );

                // Stop and exit the thread if the simulation is not running anymore
                if !simulation_running.load(Relaxed) {
                    return;
                }
            });
        }

        self.simulation_running.store(true, SeqCst);
    }

    /// Stops the simulation.
    fn stop_simulation(&mut self) {
        self.simulation_running.store(false, SeqCst);
        self.end_instant = Some(instant::Instant::now());
    }

    /// Resets the spectrum by setting all channels to 0.
    fn reset_spectrum(&mut self) {
        self.channels.iter().for_each(|ch| ch.store(0, SeqCst));
        self.simulation_statistics = Arc::new(SimulationStatistics::default());
    }

    /// Adds the settings panel to the UI.
    fn ui_add_settings(&mut self, ui: &mut egui::Ui) {
        let simulation_running = self.simulation_running.load(Relaxed);
        ui.heading(RichText::from("Detector size:").size(18.0));
        egui::Grid::new("detector_size_settings_grid").show(ui, |ui| {
            ui.label("Radius: ");
            ui.add_enabled(
                !simulation_running,
                egui::Slider::new(&mut self.arguments.radius, 0.5_f64..=12.0).text("cm"),
            );
            ui.end_row();
            ui.label("Height: ");
            ui.add_enabled(
                !simulation_running,
                egui::Slider::new(&mut self.arguments.height, 0.5_f64..=12.0).text("cm"),
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
    }

    /// Adds the start and stop button to the UI, depending on the simulation state.
    fn ui_add_start_stop_button(&mut self, ui: &mut egui::Ui) {
        let simulation_running = self.simulation_running.load(Relaxed);

        ui.vertical_centered(|ui| {
            ui.add_space(15.0);

            let start_stop_text = if simulation_running {
                "Stop simulation"
            } else {
                "Start simulation"
            };
            let start_button =
                egui::widgets::Button::new(RichText::new(start_stop_text).strong().size(18.0));
            if ui.add(start_button).clicked() {
                self.start_stop_simulation();
            }
            if simulation_running {
                ui.add_space(5.0);
                ui.add(egui::Spinner::default());
            }
        });
    }

    /// Adds a panel to the UI, displaying the simulation statistics if there are any.
    fn ui_add_simulation_info(&mut self, ui: &mut egui::Ui) {
        let hits = self.channels.iter().map(|ch| ch.load(Relaxed)).sum::<u64>();
        let simulation_running = self.simulation_running.load(Relaxed);

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
        // Load simulation statistics from the atomic variables
        let total_photon_count = self.simulation_statistics.total_photon_count.load(Relaxed);
        let photon_count_with_detector = self
            .simulation_statistics
            .photon_count_with_detector
            .load(Relaxed);
        let total_efficiency_sum = self
            .simulation_statistics
            .total_efficiency_sum
            .load(Relaxed);
        let total_efficiency_square_sum = self
            .simulation_statistics
            .total_efficiency_square_sum
            .load(Relaxed);
        let total_efficiency = total_efficiency_sum / total_photon_count as f64;
        let total_efficiency_deviation = (total_efficiency_square_sum / total_photon_count as f64
            - total_efficiency * total_efficiency)
            .sqrt();
        let intrinsic_efficiency_sum = self
            .simulation_statistics
            .intrinsic_efficiency_sum
            .load(Relaxed);
        let intrinsic_efficiency_square_sum = self
            .simulation_statistics
            .intrinsic_efficiency_square_sum
            .load(Relaxed);
        let intrinsic_efficiency = intrinsic_efficiency_sum / photon_count_with_detector as f64;
        let intrinsic_efficiency_deviation = (intrinsic_efficiency_square_sum
            / photon_count_with_detector as f64
            - intrinsic_efficiency * intrinsic_efficiency)
            .sqrt();
        egui::Grid::new("simulation_info").show(ui, |ui| {
            ui.add(egui::widgets::Label::new(
                RichText::new("").strong().size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("All photons emitted:").strong().size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(total_photon_count.to_string())
                    .strong()
                    .size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("Photons into detector:").strong().size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(photon_count_with_detector.to_string())
                    .strong()
                    .size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("Hits:").strong().size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(hits.to_string()).strong().size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("Hitrate:").strong().size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(rate_text).strong().size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("").strong().size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("Total efficiency:").strong().size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(format!("{:.5}%", total_efficiency * 100.0))
                    .strong()
                    .size(18.0),
            ));
            ui.end_row();
            ui.add(egui::widgets::Label::new(
                RichText::new("\tdeviation").size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(format!("{:.5}%", total_efficiency_deviation * 100.0))
                    .strong()
                    .size(18.0),
            ));
            ui.end_row();

            ui.add(egui::widgets::Label::new(
                RichText::new("Intrinsic efficiency:").strong().size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(format!("{:.5}%", intrinsic_efficiency * 100.0))
                    .strong()
                    .size(18.0),
            ));
            ui.end_row();
            ui.add(egui::widgets::Label::new(
                RichText::new("\tdeviation").size(18.0),
            ));
            ui.add(egui::widgets::Label::new(
                RichText::new(format!("{:.5}%", intrinsic_efficiency_deviation * 100.0))
                    .strong()
                    .size(18.0),
            ));
            ui.end_row();
        });

        // Here is the code for exporting the simulation data to a file
        #[cfg(not(target_arch = "wasm32"))]
        if !simulation_running {
            ui.vertical_centered(|ui| {
                let export_button =
                    egui::widgets::Button::new(RichText::new("Export to JSON").strong().size(18.0));
                ui.add_space(15.0);
                if ui.add(export_button).clicked() {
                    let channel_width = self.get_max_energy() / (self.channels.len() as f64);
                    // Make data into a JSON object
                    let export_data = json::object! {
                        counts: self.channels.iter().map(|channel| {
                            channel.load(Relaxed)
                        }).collect::<Vec<_>>(),
                        energies: (0..self.channels.len()).map(|i| {
                            i as f64 * channel_width + channel_width / 2.0
                        }).collect::<Vec<_>>(),
                        channel_width: channel_width,
                        parameters: {
                            fwhm: self.arguments.fwhm,
                            emitter_position: {
                                x: self.arguments.rx,
                                y: self.arguments.ry,
                                z: self.arguments.rz,
                            },
                            photon_energy: self.arguments.energy,
                            detector_density: self.arguments.density,
                            detector_radius: self.arguments.radius,
                            detector_height: self.arguments.height,
                        },
                        statistics: {
                            total_photon_count:total_photon_count,
                            photon_count_with_detector:photon_count_with_detector,
                            hits:hits,
                            total_efficiency:total_efficiency,
                            total_efficiency_deviation:total_efficiency_deviation,
                            intrinsic_efficiency:intrinsic_efficiency,
                            intrinsic_efficiency_deviation:intrinsic_efficiency_deviation,
                        }
                    }
                    .to_string();

                    // Open file dialog to select a file to save to
                    let path = rfd::FileDialog::new()
                        .add_filter("json", &["json"])
                        .set_file_name("result.json")
                        .save_file();

                    // Write the data to the file
                    if let Some(path) = path {
                        let mut file = std::fs::File::create(path).unwrap();
                        file.write_all(export_data.as_bytes())
                            .expect("Writing the file to disk failed");
                    };
                };
            });
        }
    }

    /// This function adds the plot of the spectrum to the UI.
    fn plot_spectrum(&mut self, ui: &mut egui::Ui) {
        let channel_width = self.get_max_energy() / self.channel_number as f64;

        // Different plots based on if we want logarithmic or linear scale
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
    }

    /// This function executes the simulation on the WASM build.
    /// It is called from the main loop.
    /// It tries to maintain a constant framerate by limiting the number of
    /// loops it does based on the measured time.
    #[cfg(target_arch = "wasm32")]
    fn run_simulation_wasm(&self) {
        let photon_emit = PhotonEmitter::from_params(
            self.arguments.radius,
            self.arguments.height,
            self.arguments.rx,
            self.arguments.ry,
            self.arguments.rz,
        );
        let start_time = instant::Instant::now();
        let max_energy = self.get_max_energy();
        let deviation = self.get_standard_deviation_from_fwhm();
        let pos_vector =
            Vector::<f64>::new(self.arguments.rx, self.arguments.ry, self.arguments.rz);
        loop {
            run_simulation_cycles(
                self.arguments.energy,
                pos_vector,
                &photon_emit,
                deviation,
                max_energy,
                &self.channels,
                &self.simulation_statistics,
            );
            let end_time = instant::Instant::now();
            // If it took too much time (more than 30 ms), then break the loop
            if (end_time - start_time).as_micros() > 30_000 {
                break;
            }
        }
    }
}

impl eframe::App for MyApp {
    /// This function updates the GUI and runs the simulation in the WASM build.
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        ctx.set_visuals(egui::style::Visuals::dark());
        egui::SidePanel::right("main_settings_panel")
            .resizable(false)
            .show(ctx, |ui| {
                // Render the settings panel
                ui.heading(
                    RichText::from("Particle transport simulation")
                        .size(24.0)
                        .strong(),
                );

                self.ui_add_settings(ui);
                self.ui_add_start_stop_button(ui);

                if self.start_instant.is_some() {
                    self.ui_add_simulation_info(ui);
                }
            });

        // Plot the spectrum
        egui::CentralPanel::default().show(ctx, |ui| {
            self.plot_spectrum(ui);
        });

        // Run the simulation in the WASM build
        #[cfg(target_arch = "wasm32")]
        if self.simulation_running.load(Relaxed) {
            self.run_simulation_wasm();
        }
    }
}

impl Default for MyApp {
    /// This function creates a new instance of the app with the default settings.
    fn default() -> Self {
        let default_channel_number = 1024;
        Self {
            arguments: SimulationArgs::default(),
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
            simulation_statistics: Arc::new(SimulationStatistics::default()),
        }
    }
}
