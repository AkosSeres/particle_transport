#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

#[cfg(not(target_arch = "wasm32"))]
fn main() {
    let mut options = eframe::NativeOptions::default();
    options.vsync = true;
    options.transparent = false;
    options.maximized = true;
    eframe::run_native(
        "Photon transport",
        options,
        Box::new(|_cc| Box::new(particle_transport::app::MyApp::default())),
    );
}
