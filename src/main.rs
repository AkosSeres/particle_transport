#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

/// This is the main function, which is called if
/// the build is not running on the web. In the web
/// build, the `start` function is called instead
/// from the `lib.rs` file.
#[cfg(not(target_arch = "wasm32"))]
fn main() {
    // Create app instance and initialize
    let app_instance = particle_transport::app::MyApp::default();
    app_instance.init();

    // Set `egui`/`eframe` options
    let mut options = eframe::NativeOptions::default();
    options.vsync = true;
    options.transparent = false;
    options.maximized = true;

    // Run app
    eframe::run_native(
        "Photon transport",
        options,
        Box::new(|_cc| Box::new(app_instance)),
    );
}
