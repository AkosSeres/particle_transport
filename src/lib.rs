#![warn(clippy::all, rust_2018_idioms)]

pub mod app;
pub mod photon;
pub mod photon_emitter;
pub mod rand_gen;
pub mod vec3;

pub use eframe;

// ----------------------------------------------------------------------------
// When compiling for web:
#[cfg(target_arch = "wasm32")]
use eframe::wasm_bindgen::{self, prelude::*};

/// This is the entry-point for all the web-assembly.
/// This is called once from the HTML.
/// It loads the app, installs some callbacks, then returns.
/// You can add more callbacks like this if you want to call in to your code.
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn start(canvas_id: &str) -> Result<(), eframe::wasm_bindgen::JsValue> {
    // Make sure panics are logged using `console.error`.
    console_error_panic_hook::set_once();

    // Redirect tracing to console.log and friends:
    tracing_wasm::set_as_global_default();

    // Create the app instance
    let app_instance = app::MyApp::default();

    // Set default detector values
    app_instance.init();

    // Set the seed for random number generation
    use crate::rand_gen::set_rng_seed;
    set_rng_seed(instant::now() as u64);

    // Start web app
    eframe::start_web(canvas_id, Box::new(|cc| Box::new(app_instance)))
}
